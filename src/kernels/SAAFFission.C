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
  params.addRequiredCoupledVar(
      "group_flux_ordinates",
      "The angular flux ordinates for all groups. These variables are used to inform the coupleable "
      "interface (to setup Jacobians), and aren't actually required for any calculations.");
  params.addRequiredCoupledVar(
      "group_scalar_fluxes",
      "The scalar fluxes (zero'th moments of the angular fluxes) for all spectral energy groups. We use "
      "these moments to compute the fission term");
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

  const unsigned int num_ord = coupledComponents("group_flux_ordinates");
  if (num_ord != _aq.totalOrder() * _num_groups)
    mooseError("Mismatch between the angular flux ordinates and quadrature set.");

  // Fetch the flux ordinate derivatives.
  for (unsigned int i = 0; i < num_ord; ++i)
  {
    unsigned int g = i / _aq.totalOrder();
    unsigned int n = i - g * _aq.totalOrder();
    _jvar_map.emplace(coupled("group_flux_ordinates", i), std::make_pair(g, n));
  }

  const unsigned int num_scal = coupledComponents("group_scalar_fluxes");
  if (num_scal != _num_groups)
    mooseError("Mismatch between the number of scalar fluxes and the number of groups.");

  _group_scalar_fluxes.reserve(num_scal);
  for (unsigned int i = 0u; i < num_scal; ++i)
    _group_scalar_fluxes.emplace_back(&coupledValue("group_scalar_fluxes", i));
}

Real
SAAFFission::computeQpResidual()
{
  // Quit early if no fission cross-sections or fission spectra are provided.
  if (_nu_sigma_f_g[_qp].size() == 0u || _chi_g[_qp].size() == 0u)
    return 0.0;

  Real res = 0.0;
  for (unsigned int g_prime = 0u; g_prime < _num_groups; ++g_prime)
    res += MetaPhysicL::raw_value(_nu_sigma_f_g[_qp][g_prime]) * (*_group_scalar_fluxes[g_prime])[_qp];

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
