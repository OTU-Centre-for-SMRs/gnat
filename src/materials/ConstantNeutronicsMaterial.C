#include "ConstantNeutronicsMaterial.h"

registerMooseObject("GnatApp", ConstantNeutronicsMaterial);

InputParameters
ConstantNeutronicsMaterial::validParams()
{
  auto params = AbsorbingNeutronicsMaterial::validParams();
  params.addClassDescription("Provides the neutron group velocity ($v_{g}$), "
                             "neutron group absorption cross-section "
                             "($\\Sigma_{a,g}$), and the scattering cross-"
                             "section moments "
                             "($\\Sigma_{s, g', g, l}$) for "
                             "transport problems. If no scattering "
                             "cross-section moments are provided, the material "
                             "initializes all moments to 0. The properties "
                             "must be listed in decreasing order by energy.");
  params.addRequiredParam<std::vector<Real>>("group_scattering",
                                             "The group-to-group scattering "
                                             "cross-section moments for all "
                                             "energy groups.");
  params.addParam<unsigned int>("anisotropy", 0u, "The scattering anisotropy of the medium.");

  return params;
}

ConstantNeutronicsMaterial::ConstantNeutronicsMaterial(const InputParameters & parameters)
  : AbsorbingNeutronicsMaterial(parameters),
    _sigma_s_g_prime_g_l(getParam<std::vector<Real>>("group_scattering")),
    _anisotropy(getParam<unsigned int>("anisotropy")),
    _max_moments((_anisotropy + 1u) * _num_groups * _num_groups)
{
  // Warn the user if more parameters have been provided than required.
  if (_sigma_s_g_prime_g_l.size() > _max_moments)
  {
    mooseWarning("More scattering moments have been provided than possibly "
                 "supported with the given maximum anisotropy and number of "
                 "groups. The vector will be truncated.");
  }

  // Error if the user did not provide enough parameters.
  if (_sigma_s_g_prime_g_l.size() < _max_moments && _sigma_s_g_prime_g_l.size() != 0u)
  {
    mooseError("Not enough scattering cross-section moments have been "
               "provided.");
  }

  // Resize the moments vector to 0.0 if no scattering cross-section moments
  // are provided. Required for cross-compatability with SourceNeutronicsMaterial.
  if (_sigma_s_g_prime_g_l.size() == 0u)
    _sigma_s_g_prime_g_l.resize(_max_moments, 0.0);

  // Compute the out-scattering cross-section. This is the sum of the 0'th
  // moments of the scattering cross-sections from the current group into all
  // other groups (excluding the current group).
  _sigma_s_g.resize(_num_groups, 0.0);
  _sigma_s_g_g.resize(_num_groups, 0.0);
  for (unsigned int g = 0u; g < _num_groups; ++g)
  {
    for (unsigned int g_prime = 0u; g_prime < _num_groups; ++g_prime)
      _sigma_s_g[g] +=
          _sigma_s_g_prime_g_l[g_prime * _num_groups * (_anisotropy + 1u) + g * (_anisotropy + 1u)];

    _sigma_s_g_g[g] =
        _sigma_s_g_prime_g_l[g * _num_groups * (_anisotropy + 1u) + g * (_anisotropy + 1u)];
  }

  // Compute the total cross-section.
  _sigma_t_g.resize(_num_groups, 0.0);
  for (unsigned int g = 0u; g < _num_groups; ++g)
    _sigma_t_g[g] = _sigma_s_g[g] + _sigma_a_g[g];

  // Recompute the neutron diffusion coefficients to account for scattering (transport
  // approximation).
  // Sum the first order Legendre out-scattering cross-sections for all groups.
  std::vector<Real> _sigma_s_g_prime_g_1;
  _sigma_s_g_prime_g_1.resize(_num_groups, 0.0);
  if (_anisotropy > 0u)
  {
    for (unsigned int g = 0u; g < _num_groups; ++g)
    {
      for (unsigned int g_prime = 0u; g_prime < _num_groups; ++g_prime)
      {
        _sigma_s_g_prime_g_1[g] += _sigma_s_g_prime_g_l[g_prime * _num_groups * (_anisotropy + 1u) +
                                                        g * (_anisotropy + 1u) + 1u];
      }
    }
  }

  if (_is_diffusion)
  {
    _diffusion_g.resize(_num_groups, 0.0);
    bool warning = false;
    for (unsigned int g = 0u; g < _diffusion_g.size(); ++g)
    {
      if (3.0 * (_sigma_a_g[g] + _sigma_s_g[g] - _sigma_s_g_prime_g_1[g]) <
          3.0 * libMesh::TOLERANCE)
      {
        _diffusion_g[g] = 1.0 / (3.0 * libMesh::TOLERANCE);
        warning = true;
      }
      else
        _diffusion_g[g] = 1.0 / (3.0 * (_sigma_a_g[g] + _sigma_s_g[g] - _sigma_s_g_prime_g_1[g]));
    }

    if (warning)
      mooseWarning("3.0 * (_sigma_a_g[g] + _sigma_s_g[g] - _sigma_s_g_prime_g_1[g]) < "
                   "3.0 * libMesh::TOLERANCE for the "
                   "provided cross-section(s). Using a diffusion coefficient of 1 / "
                   "3.0 * libMesh::TOLERANCE for those values.");
  }

#ifdef DEBUG_CROSS_SECTIONS
  // Output material properties for verification.
  for (unsigned int g = 0u; g < _num_groups; ++g)
  {
    _console << COLOR_GREEN << "Group " << g << " cross-sections for " << _name << ":\n"
             << COLOR_DEFAULT;
    _console << "_sigma_t_g              = " << _sigma_t_g[g] << "\n";
    _console << "_sigma_a_g              = " << _sigma_a_g[g] << "\n";
    _console << "_sigma_s_g              = " << _sigma_s_g[g] << "\n";
    _console << "_sigma_s_g + _sigma_a_g = " << _sigma_s_g[g] + _sigma_a_g[g] << "\n";
    _console << "_sigma_s_g_g            = " << _sigma_s_g_g[g] << "\n";
    _console << "IsConservative:           " << std::boolalpha
             << ((_sigma_s_g[g] + _sigma_a_g[g]) <= _sigma_t_g[g]) << std::endl;
  }
#endif
}

void
ConstantNeutronicsMaterial::computeQpProperties()
{
  EmptyNeutronicsMaterial::computeQpProperties();

  // SAAF tau.
  if (_is_saaf)
  {
    (*_mat_saaf_tau)[_qp].resize(_num_groups, 0.0);

    auto h = _current_elem->hmin();
    for (unsigned int g = 0; g < _num_groups; ++g)
    {
      if (_sigma_t_g[g] * _saaf_c * h >= _saaf_eta)
        (*_mat_saaf_tau)[_qp][g] = 1.0 / (_sigma_t_g[g] * _saaf_c);
      else
        (*_mat_saaf_tau)[_qp][g] = h / _saaf_eta;
    }
  }

  // Removal cross-sections.
  _mat_sigma_t_g[_qp].resize(_num_groups, 0.0);
  for (unsigned int g = 0; g < _num_groups; ++g)
  {
    // Have to sum the absorption and group g cross-section to form the
    // total cross-section.
    _mat_sigma_t_g[_qp][g] = _sigma_t_g[g];
  }

  if (_is_diffusion)
  {
    (*_mat_sigma_r_g)[_qp].resize(_num_groups, 0.0);
    (*_mat_diffusion_g)[_qp].resize(_num_groups, 0.0);
    for (unsigned int g = 0; g < _num_groups; ++g)
    {
      // Have to sum the absorption and out-scattering cross-section to form the
      // total cross-section. This sums all g -> g_prime scattering cross-sections and then
      // subtracts the g -> g cross-section.
      (*_mat_sigma_r_g)[_qp][g] = _sigma_t_g[g] - _sigma_s_g_g[g];
      (*_mat_diffusion_g)[_qp][g] = _diffusion_g[g];
    }
  }

  // Particle speeds.
  _mat_inv_v_g[_qp].resize(_num_groups, 0.0);
  if (_particle == Particletype::Neutron)
  {
    for (unsigned int g = 0; g < _num_groups; ++g)
      _mat_inv_v_g[_qp][g] = 1.0 / _v_g[g];
  }
  if (_particle == Particletype::GammaPhoton)
  {
    for (unsigned int g = 0; g < _num_groups; ++g)
      _mat_inv_v_g[_qp][g] = _inv_c_cm;
  }

  // Scattering moments and anisotropy.
  _mat_anisotropy[_qp] = _anisotropy;
  _mat_sigma_s_g_prime_g_l[_qp].resize(_max_moments, 0.0);
  for (unsigned int g = 0u; g < _max_moments; ++g)
    _mat_sigma_s_g_prime_g_l[_qp][g] = _sigma_s_g_prime_g_l[g];
}
