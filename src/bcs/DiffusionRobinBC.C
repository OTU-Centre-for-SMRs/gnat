#include "DiffusionRobinBC.h"

registerMooseObject("GnatApp", DiffusionRobinBC);

InputParameters
DiffusionRobinBC::validParams()
{
  auto params = IntegratedBC::validParams();
  params.addClassDescription(
      "Computes the Robin boundary condition for the neutron transport equation simplified with "
      "the diffusion approximation. The weak form is given as: $\\langle \\psi_{j}, "
      "\\frac{2}{e_{g}}(J_{g}^{inc} - \\frac{\\Phi}{4})\\rangle_{\\Gamma_{R}}$. This BC should not "
      "be exposed to the user, instead being enabled through a transport action.");
  params.addRequiredParam<Real>("incoming_partial_current",
                                "The incoming partial current for spectral energy group g at the "
                                "boundary. Set to zero for a vacuum boundary condition.");
  params.addParam<Real>("boundary_transport_correction",
                        1.0,
                        "The transport correction used to force the diffusion approximation to "
                        "better match a transport solution.");

  return params;
}

DiffusionRobinBC::DiffusionRobinBC(const InputParameters & parameters)
  : IntegratedBC(parameters), _e_g(0.0), _j_g_inc(getParam<Real>("incoming_partial_current"))
{
  if (getParam<Real>("boundary_transport_correction") <= 0.0)
    mooseError("The boundary transport correction cannot be equal to or less than zero.");

  // Pre-divide to reduce the number of floating point divisions per residual evaluation.
  _e_g = 1.0 / getParam<Real>("boundary_transport_correction");
}

Real
DiffusionRobinBC::computeQpResidual()
{
  return -2.0 * _test[_qp][_j] * _e_g * (_j_g_inc - (0.25 * _u[_qp]));
}

Real
DiffusionRobinBC::computeQpJacobian()
{
  return 0.5 * _test[_qp][_j] * _e_g * _phi[_j][_qp];
}
