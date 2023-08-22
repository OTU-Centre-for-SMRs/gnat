#include "SAAFNuclideDecaySource.h"

registerMooseObject("GnatApp", SAAFNuclideDecaySource);

InputParameters
SAAFNuclideDecaySource::validParams()
{
  auto params = SAAFBaseKernel::validParams();
  params.addClassDescription("Computes the source term for the SAAF "
                             "discrete ordinates transport equation, "
                             "where the isotropic source is provided by the decay of coupled "
                             "radionuclides. The weak form is given by "
                             "$-(\\psi_{j} + \\tau_{g}\\vec{\\nabla}\\psi_{j}"
                             "\\cdot\\hat{\\Omega}, \\sum_{i = 1}^{N} "
                             "\\frac{1}{4\\pi}\\lambda_{i}N_{i})$. "
                             "This kernel should not be exposed to the user, "
                             "instead being enabled through a transport action.");
  params.addRequiredCoupledVar(
      "nuclide_number_densities",
      "The number densities of the coupled radionuclides to use as the particle source.");
  params.addRequiredParam<std::vector<Real>>(
      "decay_transfer_functions",
      "The product of the decay constant of nuclide i, decay branching factor (probability of "
      "particle j being emitted from nuclide i due to decay), and the group-wise decay spectra "
      "(probability of particle j being emitted in particle energy group g).");

  return params;
}

SAAFNuclideDecaySource::SAAFNuclideDecaySource(const InputParameters & parameters)
  : SAAFBaseKernel(parameters),
    _decay_transfer_function(getParam<std::vector<Real>>("decay_transfer_functions"))
{
  const unsigned int num_coupled = coupledComponents("nuclide_number_densities");

  if (_decay_transfer_function.size() != num_coupled)
  {
    mooseError("The number of provided nuclide decay constants is not equal to the number of "
               "provided nuclide densities.");
  }

  _source_nuclide_densities.reserve(num_coupled);
  for (unsigned int i = 0; i < num_coupled; ++i)
  {
    _source_nuclide_densities.emplace_back(&coupledValue("nuclide_number_densities", i));
    _jvar_map.emplace(coupled("nuclide_number_densities", i), i);
  }
}

Real
SAAFNuclideDecaySource::computeQpResidual()
{
  Real res = 0.0;
  for (unsigned int i = 0u; i < _source_nuclide_densities.size(); ++i)
    res += (*_source_nuclide_densities[i])[_qp] * _decay_transfer_function[i];

  return -1.0 * computeQpTests() * res / (4.0 * libMesh::pi) * _symmetry_factor;
}

Real
SAAFNuclideDecaySource::computeQpOffDiagJacobian(unsigned int jvar)
{
  if (_jvar_map.count(jvar) == 0)
    return 0.0;

  auto i = _jvar_map[jvar];
  return -1.0 * computeQpTests() * _phi[_j][_qp] * _decay_transfer_function[i] /
         (4.0 * libMesh::pi) * _symmetry_factor;
}
