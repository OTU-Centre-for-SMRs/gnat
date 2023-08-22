#pragma once

#include "SAAFBaseKernel.h"

class SAAFNuclideDecaySource : public SAAFBaseKernel
{
public:
  static InputParameters validParams();

  SAAFNuclideDecaySource(const InputParameters & parameters);

protected:
  virtual Real computeQpResidual() override;
  virtual Real computeQpOffDiagJacobian(unsigned int jvar) override;

  // A vector of particle source nuclides. These are assumed to decay in a mono-energetic manner
  // isotropically. They are assumed to be provided as a number density.
  std::vector<const VariableValue *> _source_nuclide_densities;
  std::map<unsigned int, unsigned int> _jvar_map;

  // A vector of decay transfer functions.
  const std::vector<Real> & _decay_transfer_function;
}; // class ADSAAFNuclideDecaySource
