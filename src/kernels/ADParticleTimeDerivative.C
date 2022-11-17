#include "ADParticleTimeDerivative.h"

registerMooseObject("GnatApp", ADParticleTimeDerivative);

InputParameters
ADParticleTimeDerivative::validParams()
{
  auto params = ADTimeDerivative::validParams();
  params.addClassDescription(
      "Computes the time derivative term for either the "
      "discrete ordinates or diffusion approximation neutron transport "
      "equation. The weak form is given by $(\\psi_{j}, \\frac{1}{v_{g}}"
      "\\frac{\\partial}{\\partial t}\\Psi_{g, n}^{k})$ or "
      "$(\\psi_{j}, \\frac{1}{v_{g}} \\frac{\\partial}{\\partial t}"
      "\\Phi_{g, n}^{k})$. The SAAF discrete ordinate time derivative is treated separately to "
      "account for the stabilization terms. This kernel should not be exposed to the user, "
      "instead being enabled through a transport action.");
  params.addRequiredRangeCheckedParam<unsigned int>("group_index",
                                                    "group_index >= 0",
                                                    "The energy group index "
                                                    "of the current "
                                                    "angular flux.");
  return params;
}

ADParticleTimeDerivative::ADParticleTimeDerivative(const InputParameters & parameters)
  : ADTimeDerivative(parameters),
    _group_index(getParam<unsigned int>("group_index")),
    _inv_v_g(getADMaterialProperty<std::vector<Real>>("inv_v_g"))
{
}

ADReal
ADParticleTimeDerivative::precomputeQpResidual()
{
  if (_group_index >= _inv_v_g[_qp].size())
    mooseError("The group index exceeds the number of provided neutron speeds.");

  return _inv_v_g[_qp][_group_index] * ADTimeDerivative::precomputeQpResidual();
}
