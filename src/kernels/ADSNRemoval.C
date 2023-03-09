#include "ADSNRemoval.h"

registerMooseObject("GnatApp", ADSNRemoval);

InputParameters
ADSNRemoval::validParams()
{
  auto params = ADKernel::validParams();
  params.addClassDescription("Computes the collision term for the "
                             "discrete ordinates neutron transport equation. "
                             "The weak form is given by "
                             "$(\\psi_{j}, \\Sigma_{t,g} \\Psi_{g, n}^{k})$. "
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

ADSNRemoval::ADSNRemoval(const InputParameters & parameters)
  : ADKernel(parameters),
    _group_index(getParam<unsigned int>("group_index")),
    _sigma_t_g(getADMaterialProperty<std::vector<Real>>(getParam<std::string>("transport_system") +
                                                        "total_xs_g"))
{
}

ADReal
ADSNRemoval::computeQpResidual()
{
  if (_group_index >= _sigma_t_g[_qp].size())
  {
    mooseError("The group index exceeds the number of provided neutron total "
               "cross-sections.");
  }

  return _test[_i][_qp] * _sigma_t_g[_qp][_group_index] * _u[_qp];
}
