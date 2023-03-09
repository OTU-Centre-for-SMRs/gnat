#include "ADDFEMNuclideDecaySource.h"

registerMooseObject("GnatApp", ADDFEMNuclideDecaySource);

InputParameters
ADDFEMNuclideDecaySource::validParams()
{
  auto params = ADKernel::validParams();
  params.addClassDescription("Computes the source term for the "
                             "discrete ordinates transport equation, "
                             "where the isotropic source is provided by the decay of coupled "
                             "radionuclides. The weak form is given by "
                             "$-(\\psi_{j}, \\sum_{i = 1}^{N} "
                             "\\frac{1}{4\\pi}\\lambda_{i}N_{i})$. "
                             "This kernel should not be exposed to the user, "
                             "instead being enabled through a transport action.");
  params.addRequiredCoupledVar(
      "nuclide_number_densities",
      "The number densities of the coupled radionuclides to use as the particle source.");
  params.addRequiredParam<std::vector<Real>>(
      "decay_constants",
      "The decay constants for the decay reaction into the current particle group for each coupled "
      "radionuclide.");

  return params;
}

ADDFEMNuclideDecaySource::ADDFEMNuclideDecaySource(const InputParameters & parameters)
  : ADSNBaseKernel(parameters), _decay_consts(getParam<std::vector<Real>>("decay_constants"))
{
  const unsigned int num_coupled = coupledComponents("nuclide_number_densities");

  if (_decay_consts.size() != num_coupled)
  {
    mooseError("The number of provided nuclide decay constants is not equal to the number of "
               "provided nuclide densities.");
  }

  _source_nuclide_densities.reserve(num_coupled);
  for (unsigned int i = 0; i < num_coupled; ++i)
    _source_nuclide_densities.emplace_back(&adCoupledValue("nuclide_number_densities", i));
}

ADReal
ADDFEMNuclideDecaySource::computeQpResidual()
{
  ADReal res = 0.0;
  for (unsigned int i = 0u; i < _source_nuclide_densities.size(); ++i)
    res += (*_source_nuclide_densities[i])[_qp] * _decay_consts[i];

  return -1.0 * _test[_i][_qp] * res / (4.0 * M_PI) * _symmetry_factor;
}
