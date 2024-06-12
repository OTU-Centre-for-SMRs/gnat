#include "VoidTransportMaterial.h"

registerMooseObject("GnatApp", VoidTransportMaterial);

InputParameters
VoidTransportMaterial::validParams()
{
  auto params = EmptyTransportMaterial::validParams();
  params.addClassDescription("A void material for regions with no neutron "
                             "interactions. Provides the neutron group "
                             "velocity ($v_{g}$) for transport problems and "
                             "forces maximum void stabilization for the SAAF "
                             "scheme. The properties must be listed in "
                             "decreasing order by energy.");
  params.addParam<std::vector<Real>>("group_speeds",
                                     "The neutron speeds for all energy "
                                     "groups.");
  // Suppress eta and c since we force those to stabilize fully voided regions.
  params.suppressParameter<Real>("saaf_eta");
  params.suppressParameter<Real>("saaf_c");

  return params;
}

VoidTransportMaterial::VoidTransportMaterial(const InputParameters & parameters)
  : EmptyTransportMaterial(parameters)
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

  // Force to 0.5 and 1.0 to stabilize this void region.
  _saaf_eta = 0.5;
  _saaf_c = 1.0;

  // Compute neutron diffusion coefficients. TODO: Void stabilized neutron diffusion coefficients.
  if (_is_diffusion)
  {
    _diffusion_g.resize(_num_groups, 0.0);
    for (unsigned int g = 0u; g < _diffusion_g.size(); ++g)
      _diffusion_g[g] = 1.0 / (3.0 * libMesh::TOLERANCE);

    mooseWarning("The VoidTransportMaterial uses local diffusion coefficients with a value of 1 / "
                 "3 * libMesh::TOLERANCE.");
  }
}

void
VoidTransportMaterial::computeQpProperties()
{
  EmptyTransportMaterial::computeQpProperties();

  // SAAF tau.
  if (_is_saaf)
  {
    (*_mat_saaf_tau)[_qp].resize(_num_groups, 0.0);

    auto h = _current_elem->hmin();
    for (unsigned int g = 0; g < _num_groups; ++g)
      (*_mat_saaf_tau)[_qp][g] = h / _saaf_eta;
  }

  _mat_sigma_t_g[_qp].resize(_num_groups, 0.0);
  for (unsigned int i = 0; i < _num_groups; ++i)
    _mat_sigma_t_g[_qp][i] = 0.0;

  if (_is_diffusion)
  {
    (*_mat_sigma_r_g)[_qp].resize(_num_groups, 0.0);
    (*_mat_diffusion_g)[_qp].resize(_num_groups, 0.0);
    for (unsigned int i = 0; i < _num_groups; ++i)
    {
      (*_mat_sigma_r_g)[_qp][i] = 0.0;
      (*_mat_diffusion_g)[_qp][i] = _diffusion_g[i];
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
