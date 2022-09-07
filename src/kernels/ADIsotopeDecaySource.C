#include "ADIsotopeDecaySource.h"

registerMooseObject("GnatApp", ADIsotopeDecaySource);

InputParameters
ADIsotopeDecaySource::validParams()
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
  params.addRequiredParam<std::vector<Real>>("branching_factors",
                                             "Decay branching factors for all "
                                             "isotopes which decay into the "
                                             "current isotope. Must match the "
                                             "vector of coupled isotope "
                                             "densities.");
  params.addRequiredCoupledVar("isotope_densities",
                               "The isotope number densities for all isotopes "
                               "which form the current group under neutron "
                               "bombardment.");

  return params;
}

ADIsotopeDecaySource::ADIsotopeDecaySource(const InputParameters & parameters)
  : ADIsotopeBase(parameters),
    _decay_consts(getParam<std::vector<Real>>("decay_constants")),
    _branching_factors(getParam<std::vector<Real>>("branching_factors"))
{
  const unsigned int coupled_densities = coupledComponents("isotope_densities");
  if (coupled_densities != _decay_consts.size() || coupled_densities != _branching_factors.size())
  {
    mooseError("Mismatch between the number of provided decay "
               "constants/branching factors and the densities.");
  }

  _isotope_densities.reserve(coupled_densities);
  for (unsigned int i = 0u; i < coupled_densities; ++i)
    _isotope_densities.emplace_back(&adCoupledValue("isotope_densities", i));
}

ADReal
ADIsotopeDecaySource::computeQpResidual()
{
  ADReal res = 0.0;

  for (unsigned int i = 0u; i < _isotope_densities.size(); ++i)
    res += (*_isotope_densities[i])[_qp] * _branching_factors[i] * _decay_consts[i];

  return -1.0 * computeQpTests() * res;
}
