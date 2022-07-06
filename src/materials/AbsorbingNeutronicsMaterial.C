#include "AbsorbingNeutronicsMaterial.h"

registerMooseObject("GnatApp", AbsorbingNeutronicsMaterial);

InputParameters
AbsorbingNeutronicsMaterial::validParams()
{
  auto params = BaseNeutronicsMaterial::validParams();
  params.addClassDescription("Provides the neutron group velocity ($v_{g}$) "
                             "and the neutron group removal cross-section "
                             "($\\Sigma_{r,g}$) for simple transport problems. "
                             "The properties must be listed in decreasing "
                             "order by energy.");
  params.addRequiredRangeCheckedParam<unsigned int>("num_groups",
                                                    "num_groups > 0",
                                                    "The number of neutron "
                                                    "energy groups $n_{g}$ "
                                                    "the energy spectrum is "
                                                    "divided into.");
  params.addRequiredParam<std::vector<Real>>("group_speeds",
                                             "The neutron speeds $v_{g}$ for "
                                             "all energy groups. Units of $cm s^{-1}$.");
  params.addRequiredParam<std::vector<Real>>("group_removal",
                                             "The macroscopic neutron removal "
                                             "cross-sections $\\Sigma_{r,g}$)"
                                             "for all energy groups. Units of $cm^{-1}$.");

  return params;
}

AbsorbingNeutronicsMaterial::AbsorbingNeutronicsMaterial(const InputParameters & parameters)
  : BaseNeutronicsMaterial(parameters)
  , _num_groups(getParam<unsigned int>("num_groups"))
  , _v_g(getParam<std::vector<Real>>("group_speeds"))
  , _sigma_r_g(getParam<std::vector<Real>>("group_removal"))
  , _mat_v_g(declareADProperty<std::vector<Real>>("v_g"))
  , _mat_sigma_r_g(declareADProperty<std::vector<Real>>("removal_xs_g"))
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
  _mat_v_g[_qp].resize(_num_groups, 0.0);
  _mat_sigma_r_g[_qp].resize(_num_groups, 0.0);
  for (unsigned int i = 0; i < _num_groups; ++i)
  {
    _mat_v_g[_qp][i] = _v_g[i];
    _mat_sigma_r_g[_qp][i] = _sigma_r_g[i];
  }
}
