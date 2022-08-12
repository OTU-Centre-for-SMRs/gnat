#include "VoidNeutronicsMaterial.h"

registerMooseObject("GnatApp", VoidNeutronicsMaterial);

InputParameters
VoidNeutronicsMaterial::validParams()
{
  auto params = EmptyNeutronicsMaterial::validParams();
  params.addClassDescription("A void material for regions with no neutron "
                             "interactions. Provides the neutron group "
                             "velocity ($v_{g}$) for transport problems and "
                             "forces maximum void stabilization for the SAAF "
                             "scheme. The properties must be listed in "
                             "decreasing order by energy.");
  params.addRequiredParam<std::vector<Real>>("group_speeds",
                                             "The neutron speeds for all energy "
                                             "groups.");
  // Suppress eta and c since we force those to stabilize fully voided regions.
  params.suppressParameter<Real>("saaf_eta");
  params.suppressParameter<Real>("saaf_c");

  return params;
}

VoidNeutronicsMaterial::VoidNeutronicsMaterial(const InputParameters & parameters)
  : EmptyNeutronicsMaterial(parameters)
  , _v_g(getParam<std::vector<Real>>("group_speeds"))
{
  // Warn the user if they provided more properties than required.
  if (_v_g.size() > _num_groups)
  {
    mooseWarning("More neutron speeds provided than the number of groups. The "
                 "vector will be truncated.");
  }

  // Error if the user didn't provide enough properties.
  if (_v_g.size() < _num_groups)
    mooseError("Not enough neutron speeds have been provided.");

  // Force to 0.5 and 1.0 to stabilize this void region.
  _saaf_eta = 0.5;
  _saaf_c = 1.0;
}

void
VoidNeutronicsMaterial::computeQpProperties()
{
  _mat_anisotropy[_qp] = 0u;
  _mat_src_anisotropy[_qp] = 0u;

  // SAAF stabilization properties.
  _mat_saaf_eta[_qp] = _saaf_eta;
  _mat_saaf_c[_qp] = _saaf_c;

  _mat_v_g[_qp].resize(_num_groups, 0.0);
  _mat_sigma_r_g[_qp].resize(_num_groups, 0.0);
  for (unsigned int i = 0; i < _num_groups; ++i)
  {
    _mat_v_g[_qp][i] = _v_g[i];
    _mat_sigma_r_g[_qp][i] = 0.0;
  }
}
