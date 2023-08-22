#include "ADSNReflectiveBC.h"

registerMooseObject("GnatApp", ADSNReflectiveBC);

InputParameters
ADSNReflectiveBC::validParams()
{
  auto params = ADSNBaseBC::validParams();
  params.addClassDescription(
      "Computes the reflected boundary condition with a "
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
  params.addRequiredCoupledVar("psi_ref",
                               "The variable for "
                               "the reflected angular flux.");

  return params;
}

ADSNReflectiveBC::ADSNReflectiveBC(const InputParameters & parameters)
  : ADSNBaseBC(parameters), _ordinate_index(getParam<unsigned int>("ordinate_index"))
{
  if (_ordinate_index >= _aq.totalOrder())
    mooseError("The ordinates index exceeds the number of quadrature points.");

  const unsigned int num_coupled = coupledComponents("psi_ref");
  if (num_coupled != _aq.totalOrder())
    mooseError("Mismatch between the angular flux ordinates and quadrature set.");

  _flux_ordinates.reserve(num_coupled);
  for (unsigned int i = 0; i < num_coupled; ++i)
    _flux_ordinates.emplace_back(&coupledValue("psi_ref", i));
}

unsigned int
ADSNReflectiveBC::findReflectedOrdinate()
{
  auto dir = _aq.direction(_ordinate_index);
  auto refl_dir = MetaPhysicL::raw_value(dir - (2.0 * dir * _normals[_qp]) * _normals[_qp]);

  Real mu_refl = 0.0;
  Real omega_refl = 0.0;
  cartesianToSpherical(refl_dir, mu_refl, omega_refl);

  Real q_mu = 0.0;
  Real q_omega = 0.0;
  for (unsigned int i = 0u; i < _aq.totalOrder(); ++i)
  {
    cartesianToSpherical(_aq.direction(i), q_mu, q_omega);
    auto mu_diff = std::abs(q_mu - mu_refl);
    auto omega_diff = std::abs(q_omega - omega_refl);

    if (mu_diff < 0.0001 && omega_diff < 0.0001)
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
  ADReal n_dot_omega = _aq.direction(_ordinate_index) * _normals[_qp];
  if (n_dot_omega > 0.0)
    res += _u[_qp] * n_dot_omega * _test[_i][_qp];
  else
    res += u_ref[_qp] * n_dot_omega * _test[_i][_qp];

  return res;
}
