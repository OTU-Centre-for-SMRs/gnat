#include "ADIsotopeActivation.h"

registerMooseObject("GnatApp", ADIsotopeActivation);

InputParameters
ADIsotopeActivation::validParams()
{
  auto params = ADKernel::validParams();
  params.addClassDescription("Computes the neutron absorption source for the "
                             "isotope scalar transport equation: "
                             "$-( \\psi_{j}, \\sum_{i' = 1}^{I}"
                             "\\sum_{g = 1}^{G}\\sigma_{a,g,i'\\rightarrow i}"
                             "N_{i'}\\Phi_{g} )$.");
  params.addRequiredParam<unsigned int>("num_groups",
                                        "The number of spectral neutron energy "
                                        "groups.");
  params.addRequiredParam<std::vector<Real>>("group_activation",
                                             "The microscopic neutron "
                                             "activation cross-section for "
                                             "each energy group. These "
                                             "cross-sections must be listed in "
                                             "order of initial isotope first, "
                                             "then decending order in energy "
                                             "group second.");
  params.addRequiredCoupledVar("group_scalar_fluxes",
                               "The scalar neutron fluxes for each energy "
                               "group.");
  params.addRequiredCoupledVar("isotope_densities",
                               "The isotope number densities for all isotopes "
                               "which form the current group under neutron "
                               "bombardment.");

  return params;
}

ADIsotopeActivation::ADIsotopeActivation(const InputParameters & parameters)
  : ADKernel(parameters)
  , _num_groups(getParam<unsigned int>("num_groups"))
  , _sigma_a_g(getParam<std::vector<Real>>("group_activation"))
{
  const unsigned int coupled_densities = coupledComponents("isotope_densities");
  if (coupled_densities != (_sigma_a_g.size() / _num_groups))
  {
    mooseError("Mismatch between the number of provided cross-sections and the "
               "densities.");
  }

  const unsigned int coupled_fluxes = coupledComponents("group_scalar_fluxes");
  if (coupled_fluxes != (_sigma_a_g.size() / coupled_densities)
      && coupled_fluxes != _num_groups)
  {
    mooseError("Mismatch between the number of provided cross-sections and the "
               "scalar fluxes.");
  }

  _isotope_densities.reserve(coupled_densities);
  for (unsigned int i = 0u; i < coupled_densities; ++i)
    _isotope_densities.emplace_back(&adCoupledValue("isotope_densities", i));

  _group_scalar_fluxes.reserve(coupled_fluxes);
  for (unsigned int i = 0u; i < coupled_fluxes; ++i)
    _group_scalar_fluxes.emplace_back(&adCoupledValue("group_scalar_fluxes", i));
}

ADReal
ADIsotopeActivation::computeQpResidual()
{
  ADReal res = 0.0;
  ADReal iso_res = 0.0;

  // Loop over isotopes first.
  for (unsigned int i = 0u; i < _isotope_densities.size(); ++i)
  {
    // Groups second. _sigma_a_g is flat packed in memory to preserve coherency.
    for (unsigned int g = 0u; g < _num_groups; ++g)
      iso_res += (*_group_scalar_fluxes[g])[_qp] * _sigma_a_g[i * _num_groups + g];

    res += iso_res * (*_isotope_densities[i])[_qp];
    iso_res = 0.0;
  }

  return -1.0 * _test[_i][_qp] * res;
}
