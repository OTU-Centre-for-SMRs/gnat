#include "ADMassFractionNuclideAdvection.h"

registerMooseObject("GnatApp", ADMassFractionNuclideAdvection);

InputParameters
ADMassFractionNuclideAdvection::validParams()
{
  auto params = ADIsotopeBase::validParams();
  params.addClassDescription("Computes the advection term for the isotope "
                             "scalar transport equation: "
                             "$(- ( \\vec{\\nabla}\\psi_{j}\\cdot\\vec{v},"
                             "N_{i} )_{V}$. This kernel is stabilized with the "
                             "SUPG method.");

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

ADMassFractionNuclideAdvection::ADMassFractionNuclideAdvection(const InputParameters & parameters)
  : ADIsotopeBase(parameters)
{
}

ADReal
ADMassFractionNuclideAdvection::computeQpResidual()
{
  const auto qp_arg = std::make_tuple(_current_elem, _qp, _qrule);
  const auto vel = getQpVelocity();
  // SUPG stabilizing test functions.
  auto supg = _supg_tau(qp_arg, 0u) * vel * _grad_test[_i][_qp];

  // Unstabilized contribution.
  ADReal res = -1.0 * _grad_test[_i][_qp] * vel * _density(qp_arg, 0u) * _u[_qp];

  // Upwind contribution 1.
  res += supg * _density(qp_arg, 0u) * _u[_qp] * computeQpVelDivergence();

  // Upwind contribution 2.
  res +=
      supg * vel * (_density.gradient(qp_arg, 0u) * _u[_qp] + _grad_u[_qp] * _density(qp_arg, 0u));

  return res;
}

ADReal
ADMassFractionNuclideAdvection::computeQpVelDivergence()
{
  const auto qp_arg = std::make_tuple(_current_elem, _qp, _qrule);

  const auto div_x = _vel_u.gradient(qp_arg, 0u)(0u);
  const auto div_y = _vel_v ? _vel_v->gradient(qp_arg, 0u)(1u) : 0.0;
  const auto div_z = _vel_w ? _vel_w->gradient(qp_arg, 0u)(2u) : 0.0;

  return div_x + div_y + div_z;
}
