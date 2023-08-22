#include "ADMassFractionNuclideDecaySource.h"

registerMooseObject("GnatApp", ADMassFractionNuclideDecaySource);

InputParameters
ADMassFractionNuclideDecaySource::validParams()
{
  auto params = ADIsotopeBase::validParams();
  params.addClassDescription("Computes the neutron decay source for the "
                             "isotope scalar transport equation: "
                             "$-( \\psi_{j}, \\sum_{i' = 1}^{I}\\lambda_{i'}"
                             "f_{i'\\rightarrow i}N_{i'} )$.");
  params.addRequiredParam<std::vector<Real>>("decay_constants",
                                             "Radioactive decay constants for "
                                             "all isotopes which decay into "
                                             "the current isotope. Must match "
                                             "the vector of coupled isotope "
                                             "densities.");
  params.addCoupledVar("isotope_mass_fractions",
                       "The isotope mass fractions for all isotopes "
                       "which form the current group under neutron "
                       "bombardment.");

  return params;
}

ADMassFractionNuclideDecaySource::ADMassFractionNuclideDecaySource(
    const InputParameters & parameters)
  : ADIsotopeBase(parameters), _decay_consts(getParam<std::vector<Real>>("decay_constants"))
{
  const auto num_coupled_nuclides = coupledComponents("isotope_mass_fractions");
  if (num_coupled_nuclides != _decay_consts.size())
    mooseError("Mismatch between the number of provided decay "
               "constants/branching factors and the nuclide fractions. Decay constants: " +
               std::to_string(_decay_consts.size()) +
               "; nuclide fractions: " + std::to_string(num_coupled_nuclides));

  _isotope_fractions.reserve(num_coupled_nuclides);
  for (unsigned int i = 0u; i < num_coupled_nuclides; ++i)
    _isotope_fractions.emplace_back(&adCoupledValue("isotope_mass_fractions", i));
}

ADReal
ADMassFractionNuclideDecaySource::computeQpResidual()
{
  const auto qp_arg = std::make_tuple(_current_elem, _qp, _qrule);

  ADReal res = 0.0;
  for (unsigned int i = 0u; i < _isotope_fractions.size(); ++i)
    res += (*_isotope_fractions[i])[_qp] * _density(qp_arg, 0u) * _decay_consts[i];

  return -1.0 * computeQpTests() * res;
}
