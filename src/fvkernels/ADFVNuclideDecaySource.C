#include "ADFVNuclideDecaySource.h"

registerMooseObject("GnatApp", ADFVNuclideDecaySource);

InputParameters
ADFVNuclideDecaySource::validParams()
{
  auto params = FVElementalKernel::validParams();
  params.addClassDescription("Computes the radioactive decay source for the "
                             "isotope scalar transport equation: "
                             "$-( \\psi_{j}, \\sum_{i' = 1}^{I}\\lambda_{i'}"
                             "f_{i'\\rightarrow i}N_{i'} )$.");
  params.addRequiredParam<std::vector<Real>>("decay_constants",
                                             "Radioactive decay constants for "
                                             "all isotopes which decay into "
                                             "the current isotope. Must match "
                                             "the vector of coupled isotope "
                                             "densities.");
  params.addParam<std::vector<MooseFunctorName>>(
      "isotope_number_densities",
      "The isotope number densities for all isotopes "
      "which form the current group under during radioactive decay.");

  return params;
}

ADFVNuclideDecaySource::ADFVNuclideDecaySource(const InputParameters & parameters)
  : FVElementalKernel(parameters), _decay_consts(getParam<std::vector<Real>>("decay_constants"))
{
  const auto & nuclide_names = getParam<std::vector<MooseFunctorName>>("isotope_number_densities");
  if (nuclide_names.size() != _decay_consts.size())
    mooseError("Mismatch between the number of provided decay "
               "constants/branching factors and the nuclide number densities. Decay constants: " +
               std::to_string(_decay_consts.size()) +
               ", nuclide number densities: " + std::to_string(nuclide_names.size()));

  _isotope_number_densities.reserve(nuclide_names.size());
  for (const auto & nuclide : nuclide_names)
    _isotope_number_densities.emplace_back(&getFunctor<ADReal>(nuclide));
}

ADReal
ADFVNuclideDecaySource::computeQpResidual()
{
  auto elem_args = makeElemArg(_current_elem);

  ADReal res = 0.0;
  for (unsigned int i = 0u; i < _isotope_number_densities.size(); ++i)
    res += ((*_isotope_number_densities[i])(elem_args, 0u)) * _decay_consts[i];

  return -1.0 * res;
}
