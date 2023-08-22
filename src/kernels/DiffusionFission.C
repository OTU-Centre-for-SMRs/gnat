#include "DiffusionFission.h"

registerMooseObject("GnatApp", DiffusionFission);

InputParameters
DiffusionFission::validParams()
{
  auto params = Kernel::validParams();
  params.addClassDescription(
      "Computes the fission source term for the radiation transport "
      "equation (specialized for neutrons) simplified with the diffusion approximation. The weak "
      "form is given by: $-(\\psi_{j}, \\chi_{g}\\sum_{g' = "
      "1}^{G}\\nu\\Sigma_{f,g}\\Phi_{g',0,0})$. This kernel should not be exposed to the user, "
      "instead being enabled through a transport action.");
  params.addRequiredCoupledVar(
      "group_scalar_fluxes",
      "The scalar fluxes (zero'th moments of the angular fluxes) for all spectral energy groups.");
  params.addRequiredRangeCheckedParam<unsigned int>("group_index",
                                                    "group_index >= 0",
                                                    "The energy group index "
                                                    "of the current angular "
                                                    "flux.");
  params.addRequiredRangeCheckedParam<unsigned int>("num_groups",
                                                    "num_groups >= 1",
                                                    "The number of spectral "
                                                    "energy groups.");

  params.addParam<std::string>(
      "transport_system",
      "",
      "Name of the transport system which will consume the provided material properties. If one is "
      "not provided the first transport system will be used.");

  return params;
}

DiffusionFission::DiffusionFission(const InputParameters & parameters)
  : Kernel(parameters),
    _group_index(getParam<unsigned int>("group_index")),
    _num_groups(getParam<unsigned int>("num_groups")),
    _nu_sigma_f_g(getADMaterialProperty<std::vector<Real>>(
        getParam<std::string>("transport_system") + "production_xs_g")),
    _chi_g(getADMaterialProperty<std::vector<Real>>(getParam<std::string>("transport_system") +
                                                    "fission_spectra_g"))
{
  if (_group_index >= _num_groups)
    mooseError("The group index exceeds the number of energy groups.");

  const unsigned int num_coupled = coupledComponents("group_scalar_fluxes");
  if (num_coupled != _num_groups)
    mooseError("Mismatch between the number of scalar fluxes and the number of energy groups.");

  _group_scalar_fluxes.reserve(num_coupled);
  for (unsigned int i = 0; i < num_coupled; ++i)
  {
    _jvar_map.emplace(coupled("group_scalar_fluxes", i), i);
    _group_scalar_fluxes.emplace_back(&coupledValue("group_scalar_fluxes", i));
  }
}

Real
DiffusionFission::computeQpResidual()
{
  // Quit early if no fission cross-sections or fission spectra are provided.
  if (_nu_sigma_f_g[_qp].size() == 0u || _chi_g[_qp].size() == 0u)
    return 0.0;

  Real res = 0.0;
  for (unsigned int g_prime = 0u; g_prime < _num_groups; ++g_prime)
    res += MetaPhysicL::raw_value(_nu_sigma_f_g[_qp][g_prime]) *
           (*(_group_scalar_fluxes[g_prime]))[_qp];

  res *= MetaPhysicL::raw_value(_chi_g[_qp][_group_index]);
  return -1.0 * res * _test[_i][_qp];
}

Real
DiffusionFission::computeQpJacobian()
{
  // Quit early if no fission cross-sections or fission spectra are provided.
  if (_nu_sigma_f_g[_qp].size() == 0u || _chi_g[_qp].size() == 0u)
    return 0.0;

  return -1.0 * _test[_i][_qp] * MetaPhysicL::raw_value(_chi_g[_qp][_group_index]) *
         MetaPhysicL::raw_value(_nu_sigma_f_g[_qp][_group_index]) * _phi[_j][_qp];
}

Real
DiffusionFission::computeQpOffDiagJacobian(unsigned int jvar)
{
  // Quit early if no fission cross-sections or fission spectra are provided.
  if (_nu_sigma_f_g[_qp].size() == 0u || _chi_g[_qp].size() == 0u)
    return 0.0;

  auto & g_prime = _jvar_map[jvar];
  if (g_prime == _group_index)
    return 0.0;

  return -1.0 * _test[_i][_qp] * MetaPhysicL::raw_value(_chi_g[_qp][_group_index]) *
         MetaPhysicL::raw_value(_nu_sigma_f_g[_qp][g_prime]) * _phi[_j][_qp];
}
