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
  params.addParam<unsigned int>("source_anisotropy", 0u,
                                "The external source anisotropy of the medium.");

  return params;
}

SourceNeutronicsMaterial::SourceNeutronicsMaterial(const InputParameters & parameters)
  : ConstantNeutronicsMaterial(parameters)
  , _source_moments(getParam<std::vector<Real>>("group_source"))
  , _source_anisotropy(getParam<unsigned int>("source_anisotropy"))
{
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
  // SAAF stabilization properties.
  _mat_saaf_eta[_qp] = _saaf_eta;
  _mat_saaf_c[_qp] = _saaf_c;

  // Speeds and removal cross-sections.
  _mat_v_g[_qp].resize(_num_groups, 0.0);
  _mat_sigma_r_g[_qp].resize(_num_groups, 0.0);
  for (unsigned int i = 0; i < _num_groups; ++i)
  {
    _mat_v_g[_qp][i] = _v_g[i];
    // Have to sum the absorption and out-scattering cross-section to form the
    // removal cross-section.
    _mat_sigma_r_g[_qp][i] = _sigma_a_g[i] + _sigma_s_out[i];
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
