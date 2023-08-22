#include "SourceNeutronicsMaterial.h"

registerMooseObject("GnatApp", SourceNeutronicsMaterial);

InputParameters
SourceNeutronicsMaterial::validParams()
{
  auto params = ConstantNeutronicsMaterial::validParams();
  params.addClassDescription("Provides the neutron group velocity ($v_{g}$), "
                             "neutron group absorption cross-section "
                             "($\\Sigma_{a,g}$), the scattering cross-section "
                             "moments ($\\Sigma_{s, g', g, l}$), and the group "
                             "neutron source moments ($S_{ext,g,l,m}$) for "
                             "transport problems. If no scattering "
                             "cross-section moments are provided, the material "
                             "initializes all moments to 0. The properties "
                             "must be listed in decreasing order by energy.");
  params.addRequiredParam<std::vector<Real>>("group_source",
                                             "The external source moments for "
                                             "all energy groups.");
  params.addParam<unsigned int>(
      "source_anisotropy", 0u, "The external source anisotropy of the medium.");

  return params;
}

SourceNeutronicsMaterial::SourceNeutronicsMaterial(const InputParameters & parameters)
  : ConstantNeutronicsMaterial(parameters),
    _source_moments(getParam<std::vector<Real>>("group_source")),
    _source_anisotropy(getParam<unsigned int>("source_anisotropy"))
{
  mooseDeprecated("SourceNeutronicsMaterial is deprecated and will be removed in future versions.");

  switch (_mesh.dimension())
  {
    case 1u:
      _max_source_moments = (_source_anisotropy + 1u);
      _max_source_moments *= _num_groups;
      break;

    case 2u:
      _max_source_moments = (_source_anisotropy + 1u) * (_source_anisotropy + 2u) / 2u;
      _max_source_moments *= _num_groups;
      break;

    case 3u:
      _max_source_moments = (_source_anisotropy + 1u) * (_source_anisotropy + 1u);
      _max_source_moments *= _num_groups;
      break;

    default:
      mooseError("Unknown mesh dimensionality.");
      break;
  }

  // Warn the user if more parameters have been provided than required.
  if (_source_moments.size() > _max_source_moments)
  {
    mooseWarning("More source moments have been provided than possibly "
                 "supported with the given maximum source anisotropy and "
                 "number of groups. The vector will be truncated.");
  }

  // Error if the user did not provide enough parameters.
  if (_source_moments.size() < _max_source_moments)
    mooseError("Not enough source moments have been provided.");
}

void
SourceNeutronicsMaterial::computeQpProperties()
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

  // Speeds and removal cross-sections.
  _mat_sigma_t_g[_qp].resize(_num_groups, 0.0);
  for (unsigned int i = 0; i < _num_groups; ++i)
  {
    // Have to sum the absorption and group g cross-section to form the
    // total cross-section.
    _mat_sigma_t_g[_qp][i] = _sigma_t_g[i];
  }

  if (_is_diffusion)
  {
    (*_mat_sigma_r_g)[_qp].resize(_num_groups, 0.0);
    (*_mat_diffusion_g)[_qp].resize(_num_groups, 0.0);
    for (unsigned int i = 0; i < _num_groups; ++i)
    {
      // Have to sum the absorption and out-scattering cross-section to form the
      // total cross-section. This sums all g -> g_prime scattering cross-sections and then
      // subtracts the g -> g cross-section.
      (*_mat_sigma_r_g)[_qp][i] = _sigma_a_g[i] + _sigma_s_g[i] - _sigma_s_g_g[i];
      (*_mat_diffusion_g)[_qp][i] = _diffusion_g[i];
    }
  }

  // Particle speeds.
  _mat_inv_v_g[_qp].resize(_num_groups, 0.0);
  if (_particle == Particletype::Neutron)
  {
    for (unsigned int i = 0; i < _num_groups; ++i)
      _mat_inv_v_g[_qp][i] = 1.0 / _v_g[i];
  }
  if (_particle == Particletype::GammaPhoton)
  {
    for (unsigned int i = 0; i < _num_groups; ++i)
      _mat_inv_v_g[_qp][i] = _inv_c_cm;
  }

  // Scattering moments and anisotropy.
  _mat_anisotropy[_qp] = _anisotropy;
  _mat_sigma_s_g_prime_g_l[_qp].resize(_max_moments, 0.0);
  for (unsigned int i = 0u; i < _max_moments; ++i)
    _mat_sigma_s_g_prime_g_l[_qp][i] = _sigma_s_g_prime_g_l[i];

  // Source moments and anisotropy.
  _mat_src_anisotropy[_qp] = _source_anisotropy;
  _mat_source_moments[_qp].resize(_max_source_moments, 0.0);
  for (unsigned int i = 0u; i < _max_source_moments; ++i)
    _mat_source_moments[_qp][i] = _source_moments[i];
}
