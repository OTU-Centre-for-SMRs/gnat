#include "DiffusionApprox.h"

registerMooseObject("GnatApp", DiffusionApprox);

InputParameters
DiffusionApprox::validParams()
{
  auto params = Kernel::validParams();
  params.addClassDescription(
      "Computes the diffusion term for the "
      "neutron transport equation simplified with the diffusion approximation. "
      "The weak form is given by "
      "$(\\vec{\\nabla}\\psi_{j}, D_{g} \\vec{\\nabla}\\Phi_{g}^{k})$. "
      "This kernel should not be exposed to the user, "
      "instead being enabled through a transport action.");
  params.addRequiredRangeCheckedParam<unsigned int>("group_index",
                                                    "group_index >= 0",
                                                    "The energy group index "
                                                    "of the current angular "
                                                    "flux.");
  params.addParam<std::string>(
      "transport_system",
      "",
      "Name of the transport system which will consume the provided material properties. If one is "
      "not provided the first transport system will be used.");

  return params;
}

DiffusionApprox::DiffusionApprox(const InputParameters & parameters)
  : Kernel(parameters),
    _group_index(getParam<unsigned int>("group_index")),
    _diffusion_g(getADMaterialProperty<std::vector<Real>>(
        getParam<std::string>("transport_system") + "diffusion_g"))
{
}

Real
DiffusionApprox::computeQpResidual()
{
  return _grad_test[_i][_qp] * MetaPhysicL::raw_value(_diffusion_g[_qp][_group_index]) *
         _grad_u[_qp];
}

Real
DiffusionApprox::computeQpJacobian()
{
  return _grad_test[_i][_qp] * MetaPhysicL::raw_value(_diffusion_g[_qp][_group_index]) *
         _grad_phi[_j][_qp];
}
