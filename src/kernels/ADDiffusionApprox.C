#include "ADDiffusionApprox.h"

registerMooseObject("GnatApp", ADDiffusionApprox);

InputParameters
ADDiffusionApprox::validParams()
{
  auto params = ADKernel::validParams();
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

  return params;
}

ADDiffusionApprox::ADDiffusionApprox(const InputParameters & parameters)
  : ADKernel(parameters),
    _group_index(getParam<unsigned int>("group_index")),
    _diffusion_g(getADMaterialProperty<std::vector<Real>>("diffusion_g"))
{
}

ADReal
ADDiffusionApprox::computeQpResidual()
{
  if (_group_index >= _diffusion_g[_qp].size())
  {
    mooseError("The group index exceeds the number of provided neutron diffusion "
               "coefficients.");
  }

  return _grad_test[_i][_qp] * _diffusion_g[_qp][_group_index] * _grad_u[_qp];
}
