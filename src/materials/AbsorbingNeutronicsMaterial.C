#include "AbsorbingNeutronicsMaterial.h"

registerMooseObject("GnatApp", AbsorbingNeutronicsMaterial);

InputParameters
AbsorbingNeutronicsMaterial::validParams()
{
  auto params = EmptyNeutronicsMaterial::validParams();
  params.addClassDescription("Provides the neutron group velocity ($v_{g}$) "
                             "and the neutron group absorption cross-section "
                             "($\\Sigma_{a,g}$) for simple transport problems. "
                             "The properties must be listed in decreasing "
                             "order by energy.");
  params.addRequiredParam<std::vector<Real>>("group_speeds",
                                             "The neutron speeds for all energy "
                                             "groups.");
  params.addRequiredParam<std::vector<Real>>("group_absorption",
                                             "The macroscopic neutron absorption "
                                             "cross-sections for all energy "
                                             "groups.");

  return params;
}

AbsorbingNeutronicsMaterial::AbsorbingNeutronicsMaterial(const InputParameters & parameters)
  : EmptyNeutronicsMaterial(parameters),
    _v_g(getParam<std::vector<Real>>("group_speeds")),
    _sigma_a_g(getParam<std::vector<Real>>("group_absorption"))
{
  // Warn the user if they provided more properties than required.
  if (_v_g.size() > _num_groups)
  {
    mooseWarning("More neutron speeds provided than the number of groups. The "
                 "vector will be truncated.");
  }

  if (_sigma_a_g.size() > _num_groups)
  {
    mooseWarning("More neutron absorption cross-sections provided than the number "
                 "of groups. The vector will be truncated.");
  }

  // Error if the user didn't provide enough properties.
  if (_v_g.size() < _num_groups)
    mooseError("Not enough neutron speeds have been provided.");
  if (_sigma_a_g.size() < _num_groups)
    mooseError("Not enough neutron absorption cross-sections have been provided.");

  // Compute the neutron diffusion coefficient.
  if (_is_diffusion)
  {
    _diffusion_g.resize(_num_groups, 0.0);
    bool warning = false;
    for (unsigned int g = 0u; g < _diffusion_g.size(); ++g)
    {
      if (3.0 * _sigma_a_g[g] < 3.0 * libMesh::TOLERANCE)
      {
        _diffusion_g[g] = 1.0 / (3.0 * libMesh::TOLERANCE);
        warning = true;
      }
      else
        _diffusion_g[g] = 1.0 / (3.0 * _sigma_a_g[g]);
    }

    if (warning)
      mooseWarning(
          "3.0 * _sigma_a_g[i] < 3.0 * libMesh::TOLERANCE for the provided cross-section(s). "
          "Using a diffusion coefficient of 1 / 3.0 * libMesh::TOLERANCE for those values.");
  }
}

void
AbsorbingNeutronicsMaterial::computeQpProperties()
{
  EmptyNeutronicsMaterial::computeQpProperties();

  // SAAF tau.
  if (_is_saaf)
  {
    (*_mat_saaf_tau)[_qp].resize(_num_groups, 0.0);

    auto h = _current_elem->hmin();
    for (unsigned int g = 0; g < _num_groups; ++g)
    {
      if (_sigma_a_g[g] * _saaf_c * h >= _saaf_eta)
        (*_mat_saaf_tau)[_qp][g] = 1.0 / (_sigma_a_g[g] * _saaf_c);
      else
        (*_mat_saaf_tau)[_qp][g] = h / _saaf_eta;
    }
  }

  _mat_sigma_t_g[_qp].resize(_num_groups, 0.0);
  for (unsigned int g = 0; g < _num_groups; ++g)
    _mat_sigma_t_g[_qp][g] = _sigma_a_g[g];

  if (_is_diffusion)
  {
    (*_mat_sigma_r_g)[_qp].resize(_num_groups, 0.0);
    (*_mat_diffusion_g)[_qp].resize(_num_groups, 0.0);
    for (unsigned int g = 0; g < _num_groups; ++g)
    {
      (*_mat_sigma_r_g)[_qp][g] = _sigma_a_g[g];
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
}
