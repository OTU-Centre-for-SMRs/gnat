#include "ADNeutronTimeDerivative.h"

registerMooseObject("GnatApp", ADNeutronTimeDerivative);

InputParameters
ADNeutronTimeDerivative::validParams()
{
  auto params = ADTimeDerivative::validParams();
  params.addClassDescription("Computes the time derivative term for the "
                             "discrete ordinates neutron transport equation. "
                             "The weak form is given by $(\\psi_{j}, "
                             "\\frac{1}{v_{g}}"
                             "\\frac{\\partial}{\\partial t}\\Psi_{g, n}^{k})$. "
                             "This kernel should not be exposed to the user, "
                             "instead being enabled through a transport action.");
  params.addRequiredRangeCheckedParam<unsigned int>("group_index",
                                                    "group_index >= 0",
                                                    "The energy group index "
                                                    "$g$ of the current "
                                                    "angular flux "
                                                    "($\\Psi_{g, n}^{k}$).");
  return params;
}

ADNeutronTimeDerivative::ADNeutronTimeDerivative(const InputParameters & parameters)
  : ADTimeDerivative(parameters)
  , _group_index(getParam<unsigned int>("group_index"))
  , _v_g(getADMaterialProperty<std::vector<Real>>("group_speeds"))
{ }

ADReal
ADNeutronTimeDerivative::precomputeQpResidual()
{
  if (_group_index >= _v_g[_qp].size())
    mooseError("The group index exceeds the number of provided neutron speeds.");

  return (1.0 / _v_g[_qp][_group_index]) * ADTimeDerivative::precomputeQpResidual();
}
