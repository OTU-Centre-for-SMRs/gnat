#include "UncollidedFluxRayKernel.h"

// Local includes
#include "RayTracingStudy.h"

// MOOSE includes
#include "AuxiliarySystem.h"

#include "RealSphericalHarmonics.h"

registerMooseObject("RayTracingApp", UncollidedFluxRayKernel);

// Static mutex definition
Threads::spin_mutex UncollidedFluxRayKernel::_add_value_mutex;

InputParameters
UncollidedFluxRayKernel::validParams()
{
  auto params = IntegralRayKernelBase::validParams();
  params.addClassDescription("A ray kernel which computes the optical depth and scalar fluxes from "
                             "a point source using ray tracing.");

  params.addRequiredParam<AuxVariableName>(
      "variable",
      "The name of the array auxvariable that this RayKernel operates on. Each component will be a "
      "group-wise scalar flux.");

  params.addRequiredRangeCheckedParam<unsigned int>(
      "num_groups", "num_groups > 0", "The number of spectral energy groups.");
  params.addParam<unsigned int>("max_anisotropy",
                                0u,
                                "The maximum degree of anisotropy to evaluate. "
                                "Defaults to 0 for isotropic scattering.");
  params.addRequiredRangeCheckedParam<unsigned int>(
      "num_group_moments",
      "num_group_moments > 0",
      "The number of sspherical harmonics moments per energy group.");
  params.addParam<std::string>("source_and_weights_name",
                               "source_and_weights",
                               "The name of the ray data which houses the source intensity "
                               "and spatial weights, pre-multiplied.");
  params.addParam<std::string>(
      "transport_system",
      "",
      "Name of the transport system which will consume the provided material properties. If one is "
      "not provided the first transport system will be used.");

  return params;
}

UncollidedFluxRayKernel::UncollidedFluxRayKernel(const InputParameters & parameters)
  : IntegralRayKernelBase(parameters),
    MooseVariableInterface<RealEigenVector>(
        this,
        _fe_problem.getAuxiliarySystem()
            .getVariable(_tid, parameters.get<AuxVariableName>("variable"))
            .isNodal(),
        "variable",
        Moose::VarKindType::VAR_AUXILIARY,
        Moose::VarFieldType::VAR_FIELD_ARRAY),
    _aux(_fe_problem.getAuxiliarySystem()),
    _var(*this->mooseVariable()),
    _target_in_element(_study.getRayDataIndex("target_source_same_element")),
    _num_groups(getParam<unsigned int>("num_groups")),
    _max_eval_anisotropy(getParam<unsigned int>("max_anisotropy")),
    _num_group_moments(getParam<unsigned int>("num_group_moments")),
    _sigma_t_g(getADMaterialProperty<std::vector<Real>>(getParam<std::string>("transport_system") +
                                                        "total_xs_g"))
{
  // We do not allow RZ/RSPHERICAL because in the context of these coord
  // systems there is no way to represent a line source - we would end up
  // with a plane/surface source or a volumetric source, respectively.
  // This is also why we do not multiply by _coord[_qp] in any of the
  // integrations that follow.
  for (const auto & subdomain_id : _mesh.meshSubdomains())
    if (_fe_problem.getCoordSystem(subdomain_id) != Moose::COORD_XYZ)
      mooseError("Not valid on coordinate systems other than XYZ");

  if (_var.feType() != FEType(CONSTANT, MONOMIAL))
    paramError("variable", "Only CONSTANT MONOMIAL variables are supported");

  addMooseVariableDependency(&variable());

  // Fetch and register data indices required for group-wise calculation of the scalar flux using
  // ray-tracing. These are the group-wise source and optical depth.
  _integral_data_indices.reserve(_num_groups);
  _source_spatial_weights.reserve(_num_groups);
  for (unsigned int g = 0u; g < _num_groups; ++g)
  {
    _integral_data_indices.emplace_back(
        _study.registerRayData(integralRayDataName() + "_" + Moose::stringify(g)));
    _source_spatial_weights.emplace_back(_study.getRayDataIndex(
        getParam<std::string>("source_and_weights_name") + "_" + Moose::stringify(g), true));
  }
}

void
UncollidedFluxRayKernel::addValue(const RealEigenVector & value)
{
  Threads::spin_mutex::scoped_lock lock(_add_value_mutex);
  _var.setNodalValue(value, 0u);
  _var.add(_aux.solution());
}

void
UncollidedFluxRayKernel::onSegment()
{
  if (currentRay()->data(_target_in_element) > 0.0)
  {
    computeUncollidedFluxSourceIsTarget();
    currentRay()->setShouldContinue(false);
  }
  else
  {
    computeSegmentOpticalDepth();

    if (currentRay()->atEnd())
      computeUncollidedFluxSourceNotTarget();
  }
}

void
UncollidedFluxRayKernel::cartesianToSpherical(const RealVectorValue & direction,
                                              Real & mu,
                                              Real & omega)
{
  mu = 0.0;
  omega = 0.0;

  mu = direction(0);

  if (direction(1) > 0.0)
  {
    if (direction(2) > 0.0)
      omega += std::atan(std::abs(direction(2)) / std::abs(direction(1)));
    else if (direction(2) < 0.0)
      omega += std::atan(std::abs(direction(2)) / std::abs(direction(1))) + libMesh::pi;
  }
  else if (direction(1) < 0.0)
  {
    if (direction(2) > 0.0)
      omega += std::atan(std::abs(direction(2)) / std::abs(direction(1))) + 3.0 * libMesh::pi / 2.0;
    else if (direction(2) < 0.0)
      omega += std::atan(std::abs(direction(2)) / std::abs(direction(1))) + libMesh::pi / 2.0;
    else
      omega += libMesh::pi;
  }
  else
  {
    if (direction(2) > 0.0)
      omega += libMesh::pi / 2.0;
    else if (direction(2) < 0.0)
      omega += 3.0 * libMesh::pi / 2.0;
  }
}

// Function to compute the segment contribution to the optical depth.
void
UncollidedFluxRayKernel::computeSegmentOpticalDepth()
{
  Real integral = 0;
  for (unsigned int g = 0u; g < _num_groups; ++g)
  {
    for (_qp = 0; _qp < _q_point.size(); ++_qp)
      integral += _JxW[_qp] * MetaPhysicL::raw_value(_sigma_t_g[_qp][g]);

    // Accumulate the optical depth into the ray.
    currentRay()->data(_integral_data_indices[g]) += integral;
    integral = 0.0;
  }
}

// There is UB in this function even with addValue(val) commented out.
// - Has to do with using threads while using processors.
// - Using threading results in occasional non-deterministic segmentation faults.
// - Incorrect results even when threading is disabled, though the results are deterministic.
// - Using quadrature point zero when accessing cross-section data seems to fix the issue.
// TODO: Figure out why threaded material property access at _qp != 0 causes segmentation faults.
// Current solution where the 0th quadrature point is used works and results in the correct answer
// for target == source cells, though it is rather hacky. This might be a deeper issue than just
// this function.
void
UncollidedFluxRayKernel::computeUncollidedFluxSourceIsTarget()
{
  const auto & ray = currentRay();

  RealEigenVector val(_num_groups * _num_group_moments);
  val.setZero();

  unsigned int index = 0u;
  Real mu = 0.0;
  Real omega = 0.0;
  for (unsigned int g = 0u; g < _num_groups; ++g)
  {
    if (_mesh.dimension() == 2u)
    {
      for (unsigned int l = 0; l <= _max_eval_anisotropy; ++l)
      {
        for (int m = 0; m <= static_cast<int>(l); ++m)
        {
          if (MetaPhysicL::raw_value(_sigma_t_g[0u][g]) < libMesh::TOLERANCE)
            val[index] = ray->distance();
          else
          {
            // (1 - e^{-\Sigma_{t} * ||r_{q+} - r_{q}||}) / \Sigma_{t}
            val[index] = 1.0 / MetaPhysicL::raw_value(_sigma_t_g[0u][g]);
            val[index] *= (1.0 - std::exp(-1.0 * MetaPhysicL::raw_value(_sigma_t_g[0u][g]) *
                                          ray->distance()));
          }

          // w_{n} * w_{q} * S(r_{q'}, \hat{\Omega}_{n})
          val[index] *= ray->data(_source_spatial_weights[g]);

          // Spherical harmonics basis functions go here.
          cartesianToSpherical(ray->direction().unit(), mu, omega);
          val[index] *= RealSphericalHarmonics::evaluate(l, m, mu, omega);

          index++;
        }
      }
    }
    else
    {
      for (unsigned int l = 0; l <= _max_eval_anisotropy; ++l)
      {
        for (int m = -1 * static_cast<int>(l); m <= static_cast<int>(l); ++m)
        {
          if (MetaPhysicL::raw_value(_sigma_t_g[0u][g]) < libMesh::TOLERANCE)
            val[index] = ray->distance();
          else
          {
            // (1 - e^{-\Sigma_{t} * ||r_{q+} - r_{q}||}) / \Sigma_{t}
            val[index] = 1.0 / MetaPhysicL::raw_value(_sigma_t_g[0u][g]);
            val[index] *= (1.0 - std::exp(-1.0 * MetaPhysicL::raw_value(_sigma_t_g[0u][g]) *
                                          ray->distance()));
          }

          // w_{n} * w_{q} * S(r_{q'}, \hat{\Omega}_{n})
          val[index] *= ray->data(_source_spatial_weights[g]);

          // Spherical harmonics basis functions go here.
          cartesianToSpherical(ray->direction().unit(), mu, omega);
          val[index] *= RealSphericalHarmonics::evaluate(l, m, mu, omega);

          index++;
        }
      }
    }
  }

  // Average flux.
  val /= _current_elem->volume();

  addValue(val);
}

// Compute the uncollided flux at the destination element.
void
UncollidedFluxRayKernel::computeUncollidedFluxSourceNotTarget()
{
  const auto & ray = currentRay();

  RealEigenVector val(_num_groups * _num_group_moments);
  val.setZero();

  unsigned int index = 0u;
  Real mu = 0.0;
  Real omega = 0.0;
  for (unsigned int g = 0u; g < _num_groups; ++g)
  {
    if (_mesh.dimension() == 2u)
    {
      for (unsigned int l = 0; l <= _max_eval_anisotropy; ++l)
      {
        for (int m = 0; m <= static_cast<int>(l); ++m)
        {
          // 1 / ||r_{q'} - r_{q}||.
          val[index] = 1.0 / std::max(ray->distance(), libMesh::TOLERANCE);

          // e^{-\tau(r_{q'}, r_{q})}
          val[index] *= std::exp(-1.0 * ray->data(_integral_data_indices[g]));

          // w_{q'} * w_{q} * S(r_{q'}, r_{q'} - r_{q} / ||r_{q'} - r_{q}||)
          val[index] *= ray->data(_source_spatial_weights[g]);

          // Spherical harmonics basis functions go here.
          cartesianToSpherical(ray->direction().unit(), mu, omega);
          val[index] *= RealSphericalHarmonics::evaluate(l, m, mu, omega);

          index++;
        }
      }
    }
    else
    {
      for (unsigned int l = 0; l <= _max_eval_anisotropy; ++l)
      {
        for (int m = -1 * static_cast<int>(l); m <= static_cast<int>(l); ++m)
        {
          // 1 / ||r_{q'} - r_{q}||^2.
          val[index] = 1.0 / std::max(ray->distance() * ray->distance(),
                                      libMesh::TOLERANCE * libMesh::TOLERANCE);

          // e^{-\tau(r_{q'}, r_{q})}
          val[index] *= std::exp(-1.0 * ray->data(_integral_data_indices[g]));

          // w_{q'} * w_{q} * S(r_{q'}, r_{q'} - r_{q} / ||r_{q'} - r_{q}||)
          val[index] *= ray->data(_source_spatial_weights[g]);

          // Spherical harmonics basis functions go here.
          cartesianToSpherical(ray->direction().unit(), mu, omega);
          val[index] *= RealSphericalHarmonics::evaluate(l, m, mu, omega);

          index++;
        }
      }
    }
  }

  // Average flux.
  val /= _current_elem->volume();

  addValue(val);
}
