#include "SAAFFission.h"

registerMooseObject("GnatApp", SAAFFission);

InputParameters
SAAFFission::validParams()
{
  auto params = SAAFBaseKernel::validParams();
  params.addClassDescription(
    "Computes the fission source term for the SAAF discrete ordinates radiation transport "
    "equation (specialized for neutrons). The weak form is given by: $-(\\psi_{j} + "
    "\\tau_{g}\\vec{\\nabla}\\psi_{j}\\cdot\\hat{\\Omega}, \\frac{\\chi_{g}}{4\\pi}\\sum_{g' = "
    "1}^{G}\\nu\\Sigma_{f,g}\\Phi_{g',0,0})$. The group scalar fluxes are computed by this kernel. "
    "This kernel should not be exposed to the user, instead being enabled through a "
    "transport action.");
  params.addRequiredCoupledVar("group_flux_ordinates",
                               "The angular flux ordinates for all groups.");
  params.addRequiredRangeCheckedParam<unsigned int>("num_groups",
                                                    "num_groups >= 1",
                                                    "The number of spectral "
                                                    "energy groups.");
  return params;
}

SAAFFission::SAAFFission(const InputParameters & parameters)
  : SAAFBaseKernel(parameters),
    _num_groups(getParam<unsigned int>("num_groups")),
    _nu_sigma_f_g(getADMaterialProperty<std::vector<Real>>(
        getParam<std::string>("transport_system") + "production_xs_g")),
    _chi_g(getADMaterialProperty<std::vector<Real>>(getParam<std::string>("transport_system") +
                                                    "fission_spectra_g"))
{
  if (_group_index >= _num_groups)
  mooseError("The group index exceeds the number of energy groups.");

  if (_ordinate_index >= _aq.totalOrder())
    mooseError("The ordinates index exceeds the number of quadrature points.");

  const unsigned int num_coupled = coupledComponents("group_flux_ordinates");
  if (num_coupled != _aq.totalOrder() * _num_groups)
    mooseError("Mismatch between the angular flux ordinates and quadrature set.");

  // Fetch the flux ordinates and their derivatives.
  _group_flux_ordinates.reserve(num_coupled);
  for (unsigned int i = 0; i < num_coupled; ++i)
  {
    unsigned int g = i / _aq.totalOrder();
    unsigned int n = i - g * _aq.totalOrder();
    _jvar_map.emplace(coupled("group_flux_ordinates", i), std::make_pair(g, n));
    _group_flux_ordinates.emplace_back(&coupledValue("group_flux_ordinates", i));
  }
}

Real
SAAFFission::computeScalarFlux(unsigned int g_prime)
{
  Real moment = 0.0;
  for (unsigned int n = 0; n < _aq.totalOrder(); ++n)
    moment += (*_group_flux_ordinates[g_prime * _aq.totalOrder() + n])[_qp] * _aq.weight(n);

  return moment;
}

Real
SAAFFission::computeQpResidual()
{
  // Quit early if no fission cross-sections or fission spectra are provided.
  if (_nu_sigma_f_g[_qp].size() == 0u || _chi_g[_qp].size() == 0u)
    return 0.0;

  Real res = 0.0;
  for (unsigned int g_prime = 0u; g_prime < _num_groups; ++g_prime)
    res += MetaPhysicL::raw_value(_nu_sigma_f_g[_qp][g_prime]) * computeScalarFlux(g_prime);

  res *= MetaPhysicL::raw_value(_chi_g[_qp][_group_index]) / (4.0 * libMesh::pi) * _symmetry_factor;
  return -1.0 * res * computeQpTests();
}

Real
SAAFFission::computeQpJacobian()
{
  // Quit early if no fission cross-sections or fission spectra are provided.
  if (_nu_sigma_f_g[_qp].size() == 0u || _chi_g[_qp].size() == 0u)
    return 0.0;

  Real jac = MetaPhysicL::raw_value(_chi_g[_qp][_group_index]) / (4.0 * libMesh::pi) *
             _symmetry_factor * MetaPhysicL::raw_value(_nu_sigma_f_g[_qp][_group_index]) *
             _aq.weight(_ordinate_index) * _phi[_j][_qp];

  return -1.0 * jac * computeQpTests();
}

// TODO: Non-linear temperature and density need to be accounted for in the off-diagonal Jacobian.
Real
SAAFFission::computeQpOffDiagJacobian(unsigned int jvar)
{
  // Quit early if no fission cross-sections or fission spectra are provided.
  if (_nu_sigma_f_g[_qp].size() == 0u || _chi_g[_qp].size() == 0u)
    return 0.0;

  // Quit early if this would result in us hitting the on-diagonal Jacobian.
  auto [g_prime, n_prime] = _jvar_map[jvar];
  if (g_prime == _group_index && n_prime == _ordinate_index)
    return 0.0;

  Real jac = MetaPhysicL::raw_value(_chi_g[_qp][_group_index]) / (4.0 * libMesh::pi) *
            _symmetry_factor * MetaPhysicL::raw_value(_nu_sigma_f_g[_qp][g_prime]) *
            _aq.weight(n_prime) * _phi[_j][_qp];

  return -1.0 * computeQpTests() * jac;
}
