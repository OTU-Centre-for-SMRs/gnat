#include "ADSNRemoval.h"

registerMooseObject("GnatApp", ADSNRemoval);

InputParameters
ADSNRemoval::validParams()
{
  auto params = ADKernel::validParams();
  params.addClassDescription("Computes the removal term for the "
                             "discrete ordinates neutron transport equation. "
                             "The weak form is given by "
                             "$(\\psi_{j}, \\Sigma_{r,g} \\Psi_{g, n}^{k})$. "
                             "This kernel should not be exposed to the user, "
                             "instead being enabled through a transport action.");
  params.addRequiredRangeCheckedParam<unsigned int>("group_index",
                                                    "group_index >= 0",
                                                    "The energy group index "
                                                    "of the current angular "
                                                    "flux.");

  return params;
}

ADSNRemoval::ADSNRemoval(const InputParameters & parameters)
  : ADKernel(parameters)
  , _group_index(getParam<unsigned int>("group_index"))
  , _sigma_r_g(getADMaterialProperty<std::vector<Real>>("removal_xs_g"))
{ }

ADReal
ADSNRemoval::computeQpResidual()
{
  if (_group_index >= _sigma_r_g[_qp].size())
  {
    mooseError("The group index exceeds the number of provided neutron removal "
               "cross-sections.");
  }

  return _test[_i][_qp] * _sigma_r_g[_qp][_group_index] * _u[_qp];
}
