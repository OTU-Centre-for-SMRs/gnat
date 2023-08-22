#include "SAAFBaseDiracKernel.h"

InputParameters
SAAFBaseDiracKernel::validParams()
{
  auto params = SNBaseDiracKernel::validParams();
  params.addClassDescription("Provides stabalization parameters for SAAF "
                             "Dirac kernels, notably: $h$, $\\tau_{g}$, and "
                             "$\\phi_{j} + \\tau_{g}\\vec{\\nabla}\\phi_{j}"
                             "\\cdot\\hat{\\Omega}$. This kernel does NOT "
                             "implement computeQpResidual() or "
                             "computeQpJacobian().");
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
  params.addParam<std::string>(
      "transport_system",
      "",
      "Name of the transport system which will consume the provided material properties. If one is "
      "not provided the first transport system will be used.");

  return params;
}

SAAFBaseDiracKernel::SAAFBaseDiracKernel(const InputParameters & parameters)
  : SNBaseDiracKernel(parameters),
    _ordinate_index(getParam<unsigned int>("ordinate_index")),
    _group_index(getParam<unsigned int>("group_index")),
    _sigma_r_g(getADMaterialProperty<std::vector<Real>>(getParam<std::string>("transport_system") +
                                                        "total_xs_g")),
    _saaf_tau(getADMaterialProperty<std::vector<Real>>(getParam<std::string>("transport_system") +
                                                       "saaf_tau"))
{
}

Real
SAAFBaseDiracKernel::computeQpTests()
{
  return _test[_i][_qp] + MetaPhysicL::raw_value(_saaf_tau[_qp][_group_index]) *
                              _grad_test[_i][_qp] * _aq.direction(_ordinate_index);
}
