#include "SAAFBaseKernel.h"

InputParameters
SAAFBaseKernel::validParams()
{
  auto params = SNBaseKernel::validParams();
  params.addClassDescription("Provides stabalization parameters for SAAF "
                             "kernels, notably: $h$, $\\tau_{g}$, and "
                             "$\\psi_{j} + \\tau_{g}\\vec{\\nabla}\\psi_{j}"
                             "\\cdot\\hat{\\Omega}$. This kernel does NOT "
                             "implement computeQpResidual().");
  params.addRequiredRangeCheckedParam<unsigned int>("group_index",
                                                    "group_index >= 0",
                                                    "The energy group index "
                                                    "of the current angular "
                                                    "flux.");
  params.addRequiredRangeCheckedParam<unsigned int>("ordinate_index",
                                                    "ordinate_index >= 0",
                                                    "The discrete ordinate index "
                                                    "$n$ of the current angular "
                                                    "flux.");

  return params;
}

SAAFBaseKernel::SAAFBaseKernel(const InputParameters & parameters)
  : SNBaseKernel(parameters),
    _ordinate_index(getParam<unsigned int>("ordinate_index")),
    _group_index(getParam<unsigned int>("group_index")),
    _sigma_t_g(getADMaterialProperty<std::vector<Real>>(getParam<std::string>("transport_system") +
                                                        "total_xs_g")),
    _saaf_tau(getADMaterialProperty<std::vector<Real>>(getParam<std::string>("transport_system") +
                                                       "saaf_tau"))
{
}

Real
SAAFBaseKernel::computeQpTests()
{
  return _test[_i][_qp] + MetaPhysicL::raw_value(_saaf_tau[_qp][_group_index]) *
                              _grad_test[_i][_qp] * _aq.direction(_ordinate_index);
}
