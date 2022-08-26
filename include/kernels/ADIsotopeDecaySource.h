#pragma once

#include "ADIsotopeBase.h"

class ADIsotopeDecaySource : public ADIsotopeBase
{
public:
  static InputParameters validParams();

  ADIsotopeDecaySource(const InputParameters & parameters);

protected:
  virtual ADReal computeQpResidual() override;

  // A vector to store all of the isotope number densities which form the
  // current isotope when decaying..
  std::vector<const ADVariableValue *> _isotope_densities;

  // A vector to store all of the decay constants.
  std::vector<Real> _decay_consts;

  // A vector to store all of the branching factors.
  std::vector<Real> _branching_factors;
};
