#include "UncollidedFluxRayKernel.h"

// Local includes
#include "RayTracingStudy.h"

// MOOSE includes
#include "AuxiliarySystem.h"

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
      "variable", "The name of the auxvariable that this RayKernel operates on");

  params.addRequiredRangeCheckedParam<unsigned int>("group_index",
                                                    "group_index >= 0",
                                                    "The energy group index "
                                                    "of the current angular "
                                                    "flux.");
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
    MooseVariableInterface<Real>(this,
                                 _fe_problem.getAuxiliarySystem()
                                     .getVariable(_tid, parameters.get<AuxVariableName>("variable"))
                                     .isNodal(),
                                 "variable",
                                 Moose::VarKindType::VAR_AUXILIARY,
                                 Moose::VarFieldType::VAR_FIELD_STANDARD),
    _aux(_fe_problem.getAuxiliarySystem()),
    _var(*this->mooseVariable()),
    _integral_data_index(_study.registerRayData(integralRayDataName())),
    _source_spatial_weights(
        _study.getRayDataIndex(getParam<std::string>("source_and_weights_name"))),
    _target_in_element(_study.getRayDataIndex("target_source_same_element")),
    _group_index(getParam<unsigned int>("group_index")),
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
}

void
UncollidedFluxRayKernel::addValue(const Real value)
{
  // TODO: this is horribly inefficient. Consider caching and adding later
  Threads::spin_mutex::scoped_lock lock(_add_value_mutex);
  _var.setNodalValue(value);
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

// Function to compute the segment contribution to the optical depth.
void
UncollidedFluxRayKernel::computeSegmentOpticalDepth()
{
  Real integral = 0;
  for (_qp = 0; _qp < _q_point.size(); ++_qp)
    integral += _JxW[_qp] * MetaPhysicL::raw_value(_sigma_t_g[_qp][_group_index]);

  // Accumulate the optical depth into the ray.
  currentRay()->data(_integral_data_index) += integral;
}

void
UncollidedFluxRayKernel::computeUncollidedFluxSourceIsTarget()
{
  const auto & ray = currentRay();

  Real val = 0.0;
  if (MetaPhysicL::raw_value(_sigma_t_g[_qp][_group_index]) < libMesh::TOLERANCE)
    val += ray->data(_source_spatial_weights) * ray->distance();
  else
  {
    val += 1.0 / MetaPhysicL::raw_value(_sigma_t_g[_qp][_group_index]);
    val *= (1.0 - std::exp(-1.0 * MetaPhysicL::raw_value(_sigma_t_g[_qp][_group_index]) *
                           ray->distance()));
  }

  // w_{n} * w_{q} * S(r_{q'}, \hat{\Omega}_{n})
  val *= ray->data(_source_spatial_weights);

  // Average flux.
  val /= _current_elem->volume();

  addValue(val);
}

// Compute the uncollided flux at the destination element.
void
UncollidedFluxRayKernel::computeUncollidedFluxSourceNotTarget()
{
  const auto & ray = currentRay();

  // 1 / ||r_{q'} - r_{q}||^2.
  Real val = 1.0 / std::max(ray->distance() * ray->distance(), libMesh::TOLERANCE);

  // w_{q'} * w_{q} * S(r_{q'}, r_{q'} - r_{q} / ||r_{q'} - r_{q}||)
  val *= ray->data(_source_spatial_weights);

  // e^{-\tau(r_{q'}, r_{q})}
  val *= std::exp(-1.0 * ray->data(_integral_data_index));

  // Average flux.
  val /= _current_elem->volume();

  addValue(val);
}
