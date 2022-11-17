#include "ADSAAFBaseKernel.h"

InputParameters
ADSAAFBaseKernel::validParams()
{
  auto params = ADSNBaseKernel::validParams();
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

ADSAAFBaseKernel::ADSAAFBaseKernel(const InputParameters & parameters)
  : ADSNBaseKernel(parameters),
    _ordinate_index(getParam<unsigned int>("ordinate_index")),
    _group_index(getParam<unsigned int>("group_index")),
    _sigma_r_g(getADMaterialProperty<std::vector<Real>>("total_xs_g")),
    _saaf_eta(getADMaterialProperty<Real>("saaf_eta")),
    _saaf_c(getADMaterialProperty<Real>("saaf_c"))
{
}

ADReal
ADSAAFBaseKernel::computeQPTau()
{
  if (_group_index >= _sigma_r_g[_qp].size())
  {
    mooseError("The group index exceeds the number of provided neutron removal "
               "cross-sections.");
  }

  auto h = _current_elem->hmin();
  ADReal tau = 0.0;
  if (MetaPhysicL::raw_value(_sigma_r_g[_qp][_group_index] * _saaf_c[_qp]) * h >= _saaf_eta[_qp])
    tau = 1.0 / (_sigma_r_g[_qp][_group_index] * _saaf_c[_qp]);
  else
    tau = h / _saaf_eta[_qp];

  return tau;
}

ADReal
ADSAAFBaseKernel::computeQPTests()
{
  if (_ordinate_index >= _quadrature_set.totalOrder())
    mooseError("The ordinates index exceeds the number of quadrature points.");

  return _test[_i][_qp] +
         computeQPTau() * _grad_test[_i][_qp] * _quadrature_set.direction(_ordinate_index);
}
