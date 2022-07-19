#include "AbsorbingNeutronicsMaterial.h"

registerMooseObject("GnatApp", AbsorbingNeutronicsMaterial);

InputParameters
AbsorbingNeutronicsMaterial::validParams()
{
  auto params = EmptyNeutronicsMaterial::validParams();
  params.addClassDescription("Provides the neutron group velocity ($v_{g}$) "
                             "and the neutron group removal cross-section "
                             "($\\Sigma_{r,g}$) for simple transport problems. "
                             "The properties must be listed in decreasing "
                             "order by energy.");
  params.addRequiredParam<std::vector<Real>>("group_speeds",
                                             "The neutron speeds for all energy "
                                             "groups.");
  params.addRequiredParam<std::vector<Real>>("group_removal",
                                             "The macroscopic neutron removal "
                                             "cross-sections for all energy "
                                             "groups.");

  return params;
}

AbsorbingNeutronicsMaterial::AbsorbingNeutronicsMaterial(const InputParameters & parameters)
  : EmptyNeutronicsMaterial(parameters)
  , _v_g(getParam<std::vector<Real>>("group_speeds"))
  , _sigma_r_g(getParam<std::vector<Real>>("group_removal"))
{
  // Warn the user if they provided more properties than required.
  if (_v_g.size() > _num_groups)
  {
    mooseWarning("More neutron speeds provided than the number of groups. The "
                 "vector will be truncated.");
  }

  if (_sigma_r_g.size() > _num_groups)
  {
    mooseWarning("More neutron removal cross-sections provided than the number "
                 "of groups. The vector will be truncated.");
  }

  // Error if the user didn't provide enough properties.
  if (_v_g.size() < _num_groups)
    mooseError("Not enough neutron speeds have been provided.");
  if (_sigma_r_g.size() < _num_groups)
    mooseError("Not enough neutron removal cross-sections have been provided.");
}

void
AbsorbingNeutronicsMaterial::computeQpProperties()
{
  // Init the properties.
  if (_mat_v_g[_qp].size() != _num_groups)
  {
    _mat_v_g[_qp].resize(_num_groups, 0.0);
    _mat_sigma_r_g[_qp].resize(_num_groups, 0.0);

    for (unsigned int i = 0; i < _num_groups; ++i)
    {
      _mat_v_g[_qp][i] = _v_g[i];
      _mat_sigma_r_g[_qp][i] = _sigma_r_g[i];
    }
  }

  _mat_anisotropy[_qp] = 0u;
  _mat_src_anisotropy[_qp] = 0u;
}
