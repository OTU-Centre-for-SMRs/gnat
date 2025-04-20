#include "PropsFromVarTransportMaterial.h"

registerMooseObject("GnatApp", PropsFromVarTransportMaterial);

InputParameters
PropsFromVarTransportMaterial::validParams()
{
  auto params = EmptyTransportMaterial::validParams();
  params.addClassDescription(
      "A transport material which takes lists of cross sections (as variables) and makes them "
      "available to Gnat's particle transport solver as material properties.");

  // The different cross sections.
  params.addRequiredCoupledVar("total_xs", "The group-wise total cross sections.");
  params.addCoupledVar("scatter_xs", "The group-wise scattering matrix.");
  params.addParam<unsigned int>("anisotropy", 0u, "The scattering anisotropy of the medium.");
  params.addCoupledVar("nu_fission_xs", "The group-wise fission production cross sections.");
  params.addCoupledVar("chi", "The group-wise fission neutron spectrum.");
  params.addCoupledVar("kappa_fission", "The group-wise fission neutron heating.");
  params.addCoupledVar("inv_vel", "The group-wise inverse velocity.");
  params.addCoupledVar("absorption_xs", "The group-wise absorption cross section.");
  params.addCoupledVar("diffusion", "The group-wise particle diffusion coefficient.");

  return params;
}

PropsFromVarTransportMaterial::PropsFromVarTransportMaterial(const InputParameters & parameters)
  : EmptyTransportMaterial(parameters),
    _anisotropy(getParam<unsigned int>("anisotropy")),
    _max_moments((_anisotropy + 1u) * _num_groups * _num_groups)
{
  // Fetch the total cross sections.
  for (unsigned int g = 0; g < coupledComponents("total_xs"); ++g)
    _sigma_t_g.emplace_back(&coupledValue("total_xs", g));

  // Fetch the scattering matrix entries.
  for (unsigned int i = 0; i < coupledComponents("scatter_xs"); ++i)
    _sigma_s_g_prime_g_l.emplace_back(&coupledValue("scatter_xs", i));

  if (_has_fission)
  {
    // Fetch the neutron production cross sections.
    for (unsigned int g = 0; g < coupledComponents("nu_fission_xs"); ++g)
      _nu_sigma_f_g.emplace_back(&coupledValue("nu_fission_xs", g));

    // Fetch the chi spectrum values.
    for (unsigned int g = 0; g < coupledComponents("chi"); ++g)
      _chi_f_g.emplace_back(&coupledValue("chi", g));

    // Fetch the heating values.
    if (_has_heating)
      for (unsigned int g = 0; g < coupledComponents("kappa_fission"); ++g)
        _heating_g.emplace_back(&coupledValue("kappa_fission", g));
  }

  // Fetch the inverse velocity values.
  if (_is_transient)
    for (unsigned int g = 0; g < coupledComponents("inv_vel"); ++g)
      _inv_v_g.emplace_back(&coupledValue("inv_vel", g));

  // Fetch the diffusion coefficient values.
  if (_is_diffusion)
    for (unsigned int g = 0; g < coupledComponents("diffusion"); ++g)
      _diffusion_g.emplace_back(&coupledValue("diffusion", g));

  // Fetch the absorption cross sections.
  if (_is_diffusion)
    for (unsigned int g = 0; g < coupledComponents("absorption_xs"); ++g)
      _sigma_a_g.emplace_back(&coupledValue("absorption_xs", g));

  // Validate the provided cross sections.
  if (_sigma_t_g.size() != _num_groups && !_is_diffusion)
    paramError("total_xs",
               "Not enough total cross sections have been provided! You've provided " +
                   Moose::stringify(_sigma_t_g.size()) + " cross sections, " +
                   Moose::stringify(_num_groups) + " are required.");

  if (_sigma_s_g_prime_g_l.size() != _max_moments && !_is_diffusion)
    paramError("scatter_xs",
               "Not enough scattering cross sections have been provided! You've provided " +
                   Moose::stringify(_sigma_s_g_prime_g_l.size()) + " cross sections, " +
                   Moose::stringify(_max_moments) + " are required.");

  if (_has_fission)
  {
    if (_nu_sigma_f_g.size() != _num_groups)
      paramError(
          "nu_fission_xs",
          "Not enough neutron production cross sections have been provided! You've provided " +
              Moose::stringify(_nu_sigma_f_g.size()) + " cross sections, " +
              Moose::stringify(_num_groups) + " are required.");

    if (_chi_f_g.size() != _num_groups)
      paramError("chi",
                 "Not enough chi spectra values have been provided! You've provided " +
                     Moose::stringify(_chi_f_g.size()) + " chi values, " +
                     Moose::stringify(_num_groups) + " are required.");

    if (_heating_g.size() != _num_groups && _has_heating)
      paramError("kappa_fission",
                 "Not enough fission heating cross sections have been provided! You've provided " +
                     Moose::stringify(_heating_g.size()) + " heating values, " +
                     Moose::stringify(_num_groups) + " are required.");
  }

  if (_inv_v_g.size() != _num_groups && _is_transient)
    paramError("inv_vel",
               "Not enough inverse velocities have been provided! You've provided " +
                   Moose::stringify(_inv_v_g.size()) + " inverse velocities, " +
                   Moose::stringify(_num_groups) + " are required.");

  if (_diffusion_g.size() != _num_groups && _is_diffusion)
    paramError("diffusion",
               "Not enough diffusion coefficients have been provided! You've provided " +
                   Moose::stringify(_diffusion_g.size()) + " diffusion coefficients, " +
                   Moose::stringify(_num_groups) + " are required.");

  if (_sigma_a_g.size() != _num_groups && _is_diffusion)
    paramError("absorption_xs",
               "Not enough absorption cross sections have been provided! You've provided " +
                   Moose::stringify(_sigma_a_g.size()) + " absorption cross sections, " +
                   Moose::stringify(_num_groups) + " are required.");
}

void
PropsFromVarTransportMaterial::computeQpProperties()
{
  EmptyTransportMaterial::computeQpProperties();

  // SAAF tau.
  if (_is_saaf)
  {
    (*_mat_saaf_tau)[_qp].resize(_num_groups, 0.0);

    auto h = _current_elem->hmin();
    for (unsigned int g = 0; g < _num_groups; ++g)
    {
      if ((*(_sigma_t_g[g]))[_qp] * _saaf_c * h >= _saaf_eta)
        (*_mat_saaf_tau)[_qp][g] = 1.0 / ((*(_sigma_t_g[g]))[_qp] * _saaf_c);
      else
        (*_mat_saaf_tau)[_qp][g] = h / _saaf_eta;
    }
  }

  // Total cross sections.
  _mat_sigma_t_g[_qp].resize(_num_groups, 0.0);
  for (unsigned int g = 0; g < _num_groups; ++g)
    _mat_sigma_t_g[_qp][g] = (*(_sigma_t_g[g]))[_qp];

  // Diffusion coefficients.
  if (_is_diffusion)
  {
    (*_mat_sigma_r_g)[_qp].resize(_num_groups, 0.0);
    (*_mat_diffusion_g)[_qp].resize(_num_groups, 0.0);
    for (unsigned int g = 0; g < _num_groups; ++g)
    {
      (*_mat_sigma_r_g)[_qp][g] = (*(_sigma_a_g[g]))[_qp];
      (*_mat_diffusion_g)[_qp][g] = (*(_diffusion_g[g]))[_qp];
    }
  }

  // Particle speeds.
  if (_is_transient)
  {
    _mat_inv_v_g[_qp].resize(_num_groups, 0.0);
    for (unsigned int g = 0; g < _num_groups; ++g)
      _mat_inv_v_g[_qp][g] = (*(_inv_v_g[g]))[_qp];
  }

  // Scattering moments and anisotropy.
  _mat_anisotropy[_qp] = _anisotropy;
  _mat_sigma_s_g_prime_g_l[_qp].resize(_max_moments, 0.0);
  for (unsigned int i = 0u; i < _max_moments; ++i)
    _mat_sigma_s_g_prime_g_l[_qp][i] = (*(_sigma_s_g_prime_g_l[i]))[_qp];

  // Fission production cross-sections and spectra.
  if (_has_fission)
  {
    (*_mat_nu_sigma_f_g)[_qp].resize(_num_groups, 0.0);
    (*_mat_chi_f_g)[_qp].resize(_num_groups, 0.0);
    for (unsigned int g = 0u; g < _num_groups; ++g)
    {
      (*_mat_nu_sigma_f_g)[_qp][g] = (*(_nu_sigma_f_g[g]))[_qp];
      (*_mat_chi_f_g)[_qp][g] = (*(_chi_f_g[g]))[_qp];
    }

    // Heating values. Can only add heating if we also have fission.
    if (_has_heating)
    {
      (*_mat_heating_g)[_qp].resize(_num_groups, 0.0);
      for (unsigned int g = 0u; g < _num_groups; ++g)
        (*_mat_heating_g)[_qp][g] = (*(_heating_g[g]))[_qp];
    }
  }
}
