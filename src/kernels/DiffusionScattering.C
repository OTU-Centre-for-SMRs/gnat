#include "DiffusionScattering.h"

registerMooseObject("GnatApp", DiffusionScattering);

InputParameters
DiffusionScattering::validParams()
{
  auto params = Kernel::validParams();
  params.addClassDescription("Computes the scattering term for the "
                             "current group of the neutron transport equation simplified with the "
                             "diffusion approximation. The weak form is given by "
                             "$-(\\psi_{j}, \\sum_{g' = 1}^{G}"
                             "\\Sigma_{s,\\, g'\\rightarrow g}\\Phi_{g'})$. This kernel "
                             "should not be exposed to the user, instead being "
                             "enabled through a transport action.");
  params.addRequiredCoupledVar(
      "group_scalar_fluxes",
      "The scalar fluxes (0th degree and moment of the angular flux) for all groups.");
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

DiffusionScattering::DiffusionScattering(const InputParameters & parameters)
  : Kernel(parameters),
    _group_index(getParam<unsigned int>("group_index")),
    _num_groups(getParam<unsigned int>("num_groups")),
    _sigma_s_g_prime_g_l(getADMaterialProperty<std::vector<Real>>(
        getParam<std::string>("transport_system") + "scattering_matrix")),
    _anisotropy(getMaterialProperty<unsigned int>(getParam<std::string>("transport_system") +
                                                  "medium_anisotropy"))
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

// Compute the full scattering term.
Real
DiffusionScattering::computeQpResidual()
{
  // Quit early if no Legendre cross-section moments are provided.
  if (_sigma_s_g_prime_g_l[_qp].size() == 0u)
    return 0.0;

  Real res = 0.0;
  unsigned int scattering_index = 0u;
  for (unsigned int g_prime = 0; g_prime < _num_groups; ++g_prime)
  {
    if (g_prime == _group_index)
      continue;

    // Index into the first scattering cross-section moment.
    scattering_index =
        g_prime * _num_groups * (_anisotropy[_qp] + 1u) + _group_index * (_anisotropy[_qp] + 1u);
    res += MetaPhysicL::raw_value(_sigma_s_g_prime_g_l[_qp][scattering_index]) *
           (*_group_scalar_fluxes[g_prime])[_qp];
  }

  return -1.0 * _test[_i][_qp] * res;
}

Real
DiffusionScattering::computeQpOffDiagJacobian(unsigned int jvar)
{
  // Quit early if no Legendre cross-section moments are provided.
  if (_sigma_s_g_prime_g_l[_qp].size() == 0u)
    return 0.0;

  auto & g_prime = _jvar_map[jvar];
  if (g_prime == _group_index)
    return 0.0;

  unsigned int scattering_index =
      g_prime * _num_groups * (_anisotropy[_qp] + 1u) + _group_index * (_anisotropy[_qp] + 1u);

  return -1.0 * _test[_i][_qp] *
         MetaPhysicL::raw_value(_sigma_s_g_prime_g_l[_qp][scattering_index]) * _phi[_j][_qp];
}
