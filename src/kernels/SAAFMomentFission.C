#include "SAAFMomentFission.h"

registerMooseObject("GnatApp", SAAFMomentFission);

InputParameters
SAAFMomentFission::validParams()
{
  auto params = SAAFBaseKernel::validParams();
  params.addClassDescription(
      "Computes the fission source term for the SAAF discrete ordinates radiation transport "
      "equation (specialized for neutrons). The weak form is given by: $-(\\psi_{j} + "
      "\\tau_{g}\\vec{\\nabla}\\psi_{j}\\cdot\\hat{\\Omega}, \\frac{\\chi_{g}}{4\\pi}\\sum_{g' = "
      "1}^{G}\\nu\\Sigma_{f,g}\\Phi_{g',0,0})$. The group scalar fluxes must be provided to this "
      "kernel. This kernel should not be exposed to the user, instead being enabled through a "
      "transport action.");
  params.addRequiredCoupledVar(
      "group_scalar_fluxes",
      "The scalar fluxes (zero'th moments of the angular fluxes) for all spectral energy groups.");
  params.addRequiredRangeCheckedParam<unsigned int>("num_groups",
                                                    "num_groups >= 1",
                                                    "The number of spectral "
                                                    "energy groups.");

  return params;
}

SAAFMomentFission::SAAFMomentFission(const InputParameters & parameters)
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

  const unsigned int num_coupled = coupledComponents("group_scalar_fluxes");
  if (num_coupled != _num_groups)
    mooseError("Mismatch between the number of scalar fluxes and the number of groups.");

  _group_scalar_fluxes.reserve(num_coupled);
  for (unsigned int i = 0u; i < num_coupled; ++i)
    _group_scalar_fluxes.emplace_back(&coupledValue("group_scalar_fluxes", i));
}

Real
SAAFMomentFission::computeQpResidual()
{
  // Quit early if no fission cross-sections or fission spectra are provided.
  if (_nu_sigma_f_g[_qp].size() == 0u || _chi_g[_qp].size() == 0u)
    return 0.0;

  Real res = 0.0;
  for (unsigned int g_prime = 0u; g_prime < _num_groups; ++g_prime)
    res += MetaPhysicL::raw_value(_nu_sigma_f_g[_qp][g_prime]) *
           (*(_group_scalar_fluxes[g_prime]))[_qp];

  res *= MetaPhysicL::raw_value(_chi_g[_qp][_group_index]) / (4.0 * libMesh::pi) * _symmetry_factor;
  return -1.0 * res * computeQpTests();
}

Real
SAAFMomentFission::computeQpJacobian()
{
  // Quit early if no fission cross-sections or fission spectra are provided.
  if (_nu_sigma_f_g[_qp].size() == 0u || _chi_g[_qp].size() == 0u)
    return 0.0;

  Real jac = MetaPhysicL::raw_value(_chi_g[_qp][_group_index]) / (4.0 * libMesh::pi) *
             _symmetry_factor * MetaPhysicL::raw_value(_nu_sigma_f_g[_qp][_group_index]) *
             _aq.weight(_ordinate_index) * _phi[_j][_qp];

  return -1.0 * jac * computeQpTests();
}
