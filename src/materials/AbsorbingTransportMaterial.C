#include "AbsorbingTransportMaterial.h"

registerMooseObject("GnatApp", AbsorbingTransportMaterial);

InputParameters
AbsorbingTransportMaterial::validParams()
{
  auto params = EmptyTransportMaterial::validParams();
  params.addClassDescription("Provides the particle group velocity ($v_{g}$) "
                             "and the particle group total cross-section "
                             "($\\Sigma_{a,g}$) for simple transport problems. "
                             "The properties must be listed in decreasing "
                             "order by energy.");
  params.addRequiredParam<std::vector<Real>>("group_total",
                                             "The macroscopic particle total "
                                             "cross-sections for all energy "
                                             "groups.");
  params.addParam<std::vector<Real>>("group_speeds",
                                     "The particle speeds for all energy "
                                     "groups.");

  return params;
}

AbsorbingTransportMaterial::AbsorbingTransportMaterial(const InputParameters & parameters)
  : EmptyTransportMaterial(parameters), _sigma_t_g(getParam<std::vector<Real>>("group_total"))
{
  if (_is_transient)
  {
    if (_particle == Particletype::Neutron)
    {
      if (!isParamSetByUser("group_speeds"))
        paramError("group_speeds",
                   "The particle speeds must be provided for a transient simulation.");

      for (const auto & vel : getParam<std::vector<Real>>("group_speeds"))
        _inv_v_g.emplace_back(1.0 / vel);

      // Warn the user if they provided more properties than required.
      if (_inv_v_g.size() > _num_groups)
        mooseWarning("More neutron speeds provided than the number of groups. The "
                     "vector will be truncated.");

      // Error if the user didn't provide enough properties.
      if (_inv_v_g.size() < _num_groups)
        mooseError("Not enough neutron speeds have been provided.");
    }
    else
      _inv_v_g.resize(_num_groups, _inv_c_cm);
  }

  if (_sigma_t_g.size() > _num_groups)
    mooseWarning("More neutron absorption cross-sections provided than the number "
                 "of groups. The vector will be truncated.");

  if (_sigma_t_g.size() < _num_groups)
    mooseError("Not enough neutron absorption cross-sections have been provided.");

  // Compute the neutron diffusion coefficient.
  if (_is_diffusion)
  {
    _diffusion_g.resize(_num_groups, 0.0);
    bool warning = false;
    for (unsigned int g = 0u; g < _diffusion_g.size(); ++g)
    {
      if (3.0 * _sigma_t_g[g] < 3.0 * libMesh::TOLERANCE)
      {
        _diffusion_g[g] = 1.0 / (3.0 * libMesh::TOLERANCE);
        warning = true;
      }
      else
        _diffusion_g[g] = 1.0 / (3.0 * _sigma_t_g[g]);
    }

    if (warning)
      mooseWarning(
          "3.0 * _sigma_t_g[i] < 3.0 * libMesh::TOLERANCE for the provided cross-section(s). "
          "Using a diffusion coefficient of 1 / 3.0 * libMesh::TOLERANCE for those values.");
  }
}

void
AbsorbingTransportMaterial::computeQpProperties()
{
  EmptyTransportMaterial::computeQpProperties();

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

  _mat_sigma_t_g[_qp].resize(_num_groups, 0.0);
  for (unsigned int g = 0; g < _num_groups; ++g)
    _mat_sigma_t_g[_qp][g] = _sigma_t_g[g];

  if (_is_diffusion)
  {
    (*_mat_sigma_r_g)[_qp].resize(_num_groups, 0.0);
    (*_mat_diffusion_g)[_qp].resize(_num_groups, 0.0);
    for (unsigned int g = 0; g < _num_groups; ++g)
    {
      (*_mat_sigma_r_g)[_qp][g] = _sigma_t_g[g];
      (*_mat_diffusion_g)[_qp][g] = _diffusion_g[g];
    }
  }

  // Particle speeds.
  if (_is_transient)
  {
    _mat_inv_v_g[_qp].resize(_num_groups, 0.0);
    for (unsigned int g = 0; g < _num_groups; ++g)
      _mat_inv_v_g[_qp][g] = _inv_v_g[g];
  }
}
