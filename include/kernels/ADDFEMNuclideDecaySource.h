#pragma once

#include "ADSNBaseKernel.h"

class ADDFEMNuclideDecaySource : public ADSNBaseKernel
{
public:
  static InputParameters validParams();

  ADDFEMNuclideDecaySource(const InputParameters & parameters);

protected:
  virtual ADReal computeQpResidual() override;

  // A vector of particle source nuclides. These are assumed to decay in a mono-energetic manner
  // isotropically. They are assumed to be provided as a number density.
  std::vector<const ADVariableValue *> _source_nuclide_densities;
  // A vector of decay constants
  const std::vector<Real> & _decay_consts;
}; // class ADDFEMNuclideDecaySource
