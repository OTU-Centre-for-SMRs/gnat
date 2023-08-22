#include "DiffusionNeumannBC.h"

registerMooseObject("GnatApp", DiffusionNeumannBC);

InputParameters
DiffusionNeumannBC::validParams()
{
  auto params = IntegratedBC::validParams();
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

DiffusionNeumannBC::DiffusionNeumannBC(const InputParameters & parameters)
  : IntegratedBC(parameters), _j_g_net(getParam<Real>("net_partial_current"))
{
}

Real
DiffusionNeumannBC::computeQpResidual()
{
  return -1.0 * _test[_qp][_j] * _j_g_net;
}
