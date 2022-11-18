#include "ADDiffusionNeumannBC.h"

registerMooseObject("GnatApp", ADDiffusionNeumannBC);

InputParameters
ADDiffusionNeumannBC::validParams()
{
  auto params = ADIntegratedBC::validParams();
  params.addClassDescription(
      "Computes the Neumann boundary condition for the neutron transport equation simplified with "
      "the diffusion approximation. The weak form is given as: $\\langle \\psi_{j}, "
      "J_{g}^{net}\\rangle_{\\Gamma_{R}}$. This BC should not "
      "be exposed to the user, instead being enabled through a transport action.");
  params.addRequiredParam<Real>(
      "net_partial_current",
      "The net partial current for spectral energy group g at the "
      "boundary. Set to zero for reflective boundary conditions with no surface source.");

  return params;
}

ADDiffusionNeumannBC::ADDiffusionNeumannBC(const InputParameters & parameters)
  : ADIntegratedBC(parameters), _j_g_net(getParam<Real>("net_partial_current"))
{
}

ADReal
ADDiffusionNeumannBC::computeQpResidual()
{
  return -1.0 * _test[_qp][_j] * _j_g_net;
}
