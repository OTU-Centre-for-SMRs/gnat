#include "SUPGFunctorTestMaterial.h"

registerMooseObject("GnatApp", SUPGFunctorTestMaterial);

InputParameters
SUPGFunctorTestMaterial::validParams()
{
  auto params = FunctorMaterial::validParams();
  params += SetupInterface::validParams();
  params.set<ExecFlagEnum>("execute_on") = {EXEC_ALWAYS};
  params.addClassDescription(
      "This material provides parameters for testing isotope mass transport kernels.");

  params.addRequiredParam<NonlinearVariableName>("var_name", "The name of the variable.");
  params.addRequiredParam<Real>("diff", "The diffusion coefficient.");
  params.addRequiredParam<RealVectorValue>("vel", "The velocity.");

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

SUPGFunctorTestMaterial::SUPGFunctorTestMaterial(const InputParameters & parameters)
  : FunctorMaterial(parameters),
    _var_name(getParam<NonlinearVariableName>("var_name")),
    _diff(getParam<Real>("diff")),
    _vel(getParam<RealVectorValue>("vel"))
{
  const std::set<ExecFlagType> clearance_schedule(_execute_enum.begin(), _execute_enum.end());

  // Diffusion coefficients. This includes the turbulent contribution.
  addFunctorProperty<ADReal>(
      "isotope_diff_" + _var_name,
      [this](const auto & /*r*/, const auto & /*t*/) -> ADReal { return _diff; },
      clearance_schedule);

  // Gradient of the diffusion coefficient.
  addFunctorProperty<ADRealVectorValue>(
      "grad_isotope_diff_" + _var_name,
      [this](const auto & /*r*/, const auto & /*t*/) -> ADRealVectorValue { return 0.0; },
      clearance_schedule);

  // The stabilization parameter \tau^{SUPG}.
  addFunctorProperty<ADReal>(
      "isotope_supg_tau_" + _var_name,
      [this](const auto & /*r*/, const auto & /*t*/) -> ADReal
      {
        ADReal inv_tau_sq = 0.0;

        const auto h_min = _current_elem->hmin();
        const auto h2 = h_min * h_min;

        // A tolerance to prevent an explosion to infinity.
        const auto tol = ADRealVectorValue(libMesh::TOLERANCE * libMesh::TOLERANCE);
        const auto speed = (_vel + tol).norm();

        // Advection coontribution.
        inv_tau_sq += 4.0 * speed * speed / h2;

        // Transient contribution.
        if (_is_transient)
          inv_tau_sq += 4.0 / (_dt * _dt);

        // Diffusive contribution. There definitely has to be a better way to do this...
        inv_tau_sq += 16.0 * _diff * _diff / (h2 * h2);

        // The coefficient.
        return 1.0 / std::sqrt(inv_tau_sq);
      },
      clearance_schedule);
}
