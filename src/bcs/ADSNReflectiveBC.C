#include "ADSNReflectiveBC.h"

registerMooseObject("GnatApp", ADSNReflectiveBC);

InputParameters
ADSNReflectiveBC::validParams()
{
  auto params = ADSNBaseBC::validParams();
  params.addClassDescription("Computes the reflected boundary condition with a "
                             "weak form given by "
                             "$\\langle \\psi_{j},\\, \\hat{n}\\cdot\\hat{\\Omega}\\Psi_{r,\\, g}\\rangle_{\\Gamma_{i}}$, "
                             "$\\hat{n}\\cdot\\hat{\\Omega} \\leq 0$. "
                             "This kernel should not be exposed to the user, "
                             "instead being enabled through a transport action.");
  params.addRequiredRangeCheckedParam<unsigned int>("ordinate_index",
                                                    "ordinate_index >= 0",
                                                    "The discrete ordinate index "
                                                    "of the current angular "
                                                    "flux.");
  params.addRequiredCoupledVar("psi_ref", "The variable for "
                               "the reflected angular flux.");

  return params;
}

ADSNReflectiveBC::ADSNReflectiveBC(const InputParameters & parameters)
  : ADSNBaseBC(parameters)
  , _ordinate_index(getParam<unsigned int>("ordinate_index"))
{
  if (_ordinate_index >= _quadrature_set.totalOrder())
    mooseError("The ordinates index exceeds the number of quadrature points.");

  const unsigned int num_coupled = coupledComponents("psi_ref");
  if (num_coupled != _quadrature_set.totalOrder())
    mooseError("Mismatch between the angular flux ordinates and quadrature set.");

  _flux_ordinates.reserve(num_coupled);
  for (unsigned int i = 0; i < num_coupled; ++i)
    _flux_ordinates.emplace_back(&coupledValue("psi_ref", i));
}

unsigned int
ADSNReflectiveBC::findReflectedOrdinate()
{
  auto dir = _quadrature_set.direction(_ordinate_index);
  auto refl_dir = dir - (2.0 * dir * _normals[_qp]) * _normals[_qp];

  for (unsigned int i = 0u; i < _quadrature_set.totalOrder(); ++i)
  {
    auto x_diff = std::abs(_quadrature_set.direction(i)(0) - refl_dir(0));
    auto y_diff = std::abs(_quadrature_set.direction(i)(1) - refl_dir(1));
    auto z_diff = std::abs(_quadrature_set.direction(i)(2) - refl_dir(2));

    if (x_diff < 0.0001 && y_diff < 0.0001 && z_diff < 0.0001)
      return i;
  }

  mooseError("Reflected direction is not in the quadrature set.");
  return 0u;
}

ADReal
ADSNReflectiveBC::computeQpResidual()
{
  auto & u_ref = (*_flux_ordinates[findReflectedOrdinate()]);

  ADReal res = 0.0;
  ADReal n_dot_omega = _quadrature_set.direction(_ordinate_index) * _normals[_qp];
  if (n_dot_omega > 0.0)
    res += _u[_qp] * n_dot_omega * _test[_i][_qp];
  else
    res += u_ref[_qp] * n_dot_omega * _test[_i][_qp];

  return res;
}
