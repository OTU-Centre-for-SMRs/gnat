#include "FunctorAutoNuclideMaterial.h"

registerMooseObject("GnatApp", FunctorAutoNuclideMaterial);

InputParameters
FunctorAutoNuclideMaterial::validParams()
{
  auto params = FunctorMaterial::validParams();
  params += SetupInterface::validParams();
  params.set<ExecFlagEnum>("execute_on") = {EXEC_ALWAYS};
  params.addClassDescription(
      "This material provides parameters for isotope mass transport problems. It "
      "should not be exposed to the user, but should be automatically generated using a "
      "radionuclide system action.");

  params.addRequiredParam<NonlinearVariableName>("isotope_name",
                                                 "The name of the isotope variable.");
  params.addRequiredParam<MooseEnum>(
      "scheme",
      MooseEnum("supg_fe fv"),
      "The discretization and stabilization scheme that the nuclide system should use.");
  params.addRequiredParam<MooseEnum>("turbulence_handling",
                                     MooseEnum("none mixing-length"),
                                     "The type of diffusion coefficient to use.");

  params.addRequiredParam<Real>("radii", "The radius of the particles the field represents (cm).");
  params.addRequiredParam<MooseFunctorName>("temperature",
                                            "The temperature of the bulk fluid ($K$).");
  params.addRequiredParam<MooseFunctorName>(
      "dynamic_viscosity", "The dynamic viscosity of the bulk fluid ($g/(cm s)$).");

  params.addParam<Real>("schmidt_number",
                        0.7,
                        "The turbulent Schmidt number that relates the turbulent scalar "
                        "diffusivity to the turbulent momentum diffusivity.");
  params.addParam<MooseFunctorName>("mixing_length", "The turbulent mixing length.");

  params.addParam<MooseFunctorName>("u", "The velocity in the x direction ($cm/s$).");
  params.addParam<MooseFunctorName>("v", "The velocity in the y direction ($cm/s$).");
  params.addParam<MooseFunctorName>("w", "The velocity in the z direction ($cm/s$).");

  // Some ghosting is required for finite volume fields (neighbor information is used to compute
  // gradients).
  params.addParam<unsigned short>("ghost_layers", 2, "The number of layers of elements to ghost.");
  params.addRelationshipManager(
      "ElementSideNeighborLayers",
      Moose::RelationshipManagerType::GEOMETRIC | Moose::RelationshipManagerType::ALGEBRAIC,
      [](const InputParameters & obj_params, InputParameters & rm_params)
      {
        rm_params.set<unsigned short>("layers") = obj_params.get<unsigned short>("ghost_layers");
        rm_params.set<bool>("use_displaced_mesh") = obj_params.get<bool>("use_displaced_mesh");
      });

  return params;
}

FunctorAutoNuclideMaterial::FunctorAutoNuclideMaterial(const InputParameters & parameters)
  : FunctorMaterial(parameters),
    _mesh_dims(_subproblem.mesh().dimension()),
    _scheme(getParam<MooseEnum>("scheme").getEnum<NuclideScheme>()),
    _diff_type(getParam<MooseEnum>("turbulence_handling").getEnum<CoefficientType>()),
    _radii(getParam<Real>("radii")),
    _temperature(getFunctor<ADReal>("temperature")),
    _visc_dynamic(getFunctor<ADReal>("dynamic_viscosity")),
    _schmidt_number(getParam<Real>("schmidt_number")),
    _mixing_len(_diff_type == CoefficientType::MixingLength ? &getFunctor<ADReal>("mixing_length")
                                                            : nullptr),
    _u(_diff_type == CoefficientType::MixingLength ? &getFunctor<ADReal>("u") : nullptr),
    _v(isParamValid("v") && _diff_type == CoefficientType::MixingLength ? &getFunctor<ADReal>("v")
                                                                        : nullptr),
    _w(isParamValid("w") && _diff_type == CoefficientType::MixingLength ? &getFunctor<ADReal>("w")
                                                                        : nullptr),
    _nuc_diff(nullptr)
{
  if (_mesh_dims >= 2 && !_v && _diff_type == CoefficientType::MixingLength)
    mooseError(
        "In 2D or 3D the v component of the velocity must be supplied using the 'v' parameter.");
  if (_mesh_dims >= 3 && !_w && _diff_type == CoefficientType::MixingLength)
    mooseError("In 3D, the w component of the velocity must be supplied using the 'w' parameter.");

  const std::set<ExecFlagType> clearance_schedule(_execute_enum.begin(), _execute_enum.end());

  // Diffusion coefficients. This includes the turbulent contribution.
  addFunctorProperty<ADReal>(
      "isotope_diff_" + Moose::stringify(getParam<NonlinearVariableName>("isotope_name")),
      [this](const auto & r, const auto & t) -> ADReal
      {
        ADReal diff = (_k_boltzmann_cm * _temperature(r, t)) /
                      (6.0 * libMesh::pi * _visc_dynamic(r, t) * _radii);

        // Turbulence handling.
        if (_diff_type == CoefficientType::MixingLength)
        {
          constexpr Real offset = 1e-15; // Prevents the explosion of the derivative of sqrt(x).

          const auto grad_u = _u->gradient(r, t);
          ADReal symmetric_strain_tensor_norm = 2.0 * Utility::pow<2>(grad_u(0));
          if (_mesh_dims >= 2)
          {
            const auto grad_v = _v->gradient(r, t);
            symmetric_strain_tensor_norm +=
                2.0 * Utility::pow<2>(grad_v(1)) + Utility::pow<2>(grad_v(0) + grad_u(1));
            if (_mesh_dims >= 3)
            {
              const auto grad_w = _w->gradient(r, t);
              symmetric_strain_tensor_norm += 2.0 * Utility::pow<2>(grad_w(2)) +
                                              Utility::pow<2>(grad_u(2) + grad_w(0)) +
                                              Utility::pow<2>(grad_v(2) + grad_w(1));
            }
          }

          symmetric_strain_tensor_norm = std::sqrt(symmetric_strain_tensor_norm + offset);
          ADReal eddy_viscosity =
              symmetric_strain_tensor_norm * (*_mixing_len)(r, t) * (*_mixing_len)(r, t);

          diff += eddy_viscosity / _schmidt_number;
        }

        return diff;
      },
      clearance_schedule);

  // Declare this to make things cleaner below.
  _nuc_diff = &getFunctor<ADReal>(
      "isotope_diff_" + Moose::stringify(getParam<NonlinearVariableName>("isotope_name")));

  // The rest are only required for the SUPG finite element scheme.
  if (_scheme == NuclideScheme::SUPGFE)
  {
    // Gradient of the diffusion coefficient.
    addFunctorProperty<ADRealVectorValue>(
        "grad_isotope_diff_" + Moose::stringify(getParam<NonlinearVariableName>("isotope_name")),
        [this](const auto & r, const auto & t) -> ADRealVectorValue
        {
          ADRealVectorValue grad_diff = ADRealVectorValue(0.0);
          // Laminar contribution.
          const auto constants = _k_boltzmann_cm / (6.0 * libMesh::pi * _radii);
          const auto num = _visc_dynamic(r, t) * _temperature.gradient(r, t) -
                           _temperature(r, t) * _visc_dynamic.gradient(r, t);
          grad_diff += constants * num / (_visc_dynamic(r, t) * _visc_dynamic(r, t));

          // Turbulence handling contribution.
          // This function cannot be implemented at the moment as MOOSE does not provide 2nd
          // derivative evaluations in the Functor interface.
          // TODO: Come back once 2nd derivative computations exist.
          /*
          if (_diff_type == CoefficientType::MixingLength)
          {
          }
          */

          return grad_diff;
        },
        clearance_schedule);

    // The stabilization parameter \tau^{SUPG}.
    addFunctorProperty<ADReal>(
        "isotope_supg_tau_" + Moose::stringify(getParam<NonlinearVariableName>("isotope_name")),
        [this](const auto & r, const auto & t) -> ADReal
        {
          ADReal inv_tau_sq = 0.0;

          const auto h_min = _current_elem->hmin();
          const auto h2 = h_min * h_min;

          // A tolerance to prevent an explosion to infinity.
          const auto tol = ADRealVectorValue(libMesh::TOLERANCE * libMesh::TOLERANCE);
          const auto speed =
              (ADRealVectorValue(
                   _u ? (*_u)(r, t) : 0.0, _v ? (*_v)(r, t) : 0.0, _w ? (*_w)(r, t) : 0.0) +
               tol)
                  .norm();

          // Advection coontribution.
          inv_tau_sq += 4.0 * speed * speed / h2;

          // Transient contribution.
          if (_is_transient)
            inv_tau_sq += 4.0 / (_dt * _dt);

          // Diffusive contribution. There definitely has to be a better way to do this...
          inv_tau_sq += 16.0 * ((*_nuc_diff)(r, t)) * ((*_nuc_diff)(r, t)) / (h2 * h2);

          // The coefficient.
          return 1.0 / std::sqrt(inv_tau_sq);
        },
        clearance_schedule);
  }
}
