#include "ADDFEMTimeDerivative.h"

registerMooseObject("GnatApp", ADDFEMTimeDerivative);

InputParameters
ADDFEMTimeDerivative::validParams()
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
                                                    "of the current "
                                                    "angular flux.");
  return params;
}

ADDFEMTimeDerivative::ADDFEMTimeDerivative(const InputParameters & parameters)
  : ADTimeDerivative(parameters)
  , _group_index(getParam<unsigned int>("group_index"))
  , _v_g(getADMaterialProperty<std::vector<Real>>("v_g"))
{ }

ADReal
ADDFEMTimeDerivative::precomputeQpResidual()
{
  if (_group_index >= _v_g[_qp].size())
    mooseError("The group index exceeds the number of provided neutron speeds.");

  return (1.0 / _v_g[_qp][_group_index]) * ADTimeDerivative::precomputeQpResidual();
}
