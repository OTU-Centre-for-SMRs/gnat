#include "ADFVMassFractionNuclideActivation.h"

registerMooseObject("GnatApp", ADFVMassFractionNuclideActivation);

InputParameters
ADFVMassFractionNuclideActivation::validParams()
{
  auto params = FVElementalKernel::validParams();
  params.addClassDescription("Computes the neutron absorption source for the "
                             "isotope scalar transport equation: "
                             "$-( \\psi_{j}, \\sum_{i' = 1}^{I}"
                             "\\sum_{g = 1}^{G}\\sigma_{a,g,i'\\rightarrow i}"
                             "N_{i'}\\Phi_{g} )$.");
  params.addRequiredParam<unsigned int>("num_groups",
                                        "The number of spectral neutron energy "
                                        "groups.");
  params.addRequiredParam<MooseFunctorName>("density",
                                            "The functor name for the density of the bulk fluid.");
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
  params.addParam<std::vector<MooseFunctorName>>("isotope_mass_fractions",
                                                 "The isotope mass fractions for all isotopes "
                                                 "which form the current group under neutron "
                                                 "bombardment.");

  return params;
}

ADFVMassFractionNuclideActivation::ADFVMassFractionNuclideActivation(
    const InputParameters & parameters)
  : FVElementalKernel(parameters),
    _num_groups(getParam<unsigned int>("num_groups")),
    _density(getFunctor<ADReal>("density")),
    _sigma_a_g(getParam<std::vector<Real>>("group_activation"))
{
  const unsigned int num_coupled_fluxes = coupledComponents("group_scalar_fluxes");
  const auto & nuclide_names = getParam<std::vector<MooseFunctorName>>("isotope_mass_fractions");

  if (nuclide_names.size() != (_sigma_a_g.size() / _num_groups))
  {
    mooseError("Mismatch between the number of provided cross-sections " +
               std::to_string((_sigma_a_g.size() / _num_groups)) + " and the densities " +
               std::to_string(nuclide_names.size()) + ".");
  }

  if (num_coupled_fluxes != (_sigma_a_g.size() / nuclide_names.size()) &&
      num_coupled_fluxes != _num_groups)
  {
    mooseError("Mismatch between the number of provided cross-sections and the "
               "scalar fluxes.");
  }

  _group_scalar_fluxes.reserve(num_coupled_fluxes);
  for (unsigned int i = 0u; i < num_coupled_fluxes; ++i)
    _group_scalar_fluxes.emplace_back(&adCoupledValue("group_scalar_fluxes", i));

  _isotope_fractions.reserve(nuclide_names.size());
  for (const auto & nuclide : nuclide_names)
    _isotope_fractions.emplace_back(&getFunctor<ADReal>(nuclide));
}

ADReal
ADFVMassFractionNuclideActivation::computeQpResidual()
{
  auto elem_args = makeElemArg(_current_elem);

  ADReal res = 0.0;
  ADReal iso_res = 0.0;

  // Loop over isotopes first.
  for (unsigned int i = 0u; i < _isotope_fractions.size(); ++i)
  {
    // Groups second. _sigma_a_g is flat packed in memory to preserve coherency.
    for (unsigned int g = 0u; g < _num_groups; ++g)
      iso_res += ((*_group_scalar_fluxes[g])[_qp]) * _sigma_a_g[i * _num_groups + g];

    res += iso_res * (*_isotope_fractions[i])(elem_args, 0u);
    iso_res = 0.0;
  }

  return -1.0 * res * _density(elem_args, 0u);
}
