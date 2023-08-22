#pragma once

#include "ADIsotopeBase.h"

class ADMassFractionNuclideDecaySource : public ADIsotopeBase
{
public:
  static InputParameters validParams();

  ADMassFractionNuclideDecaySource(const InputParameters & parameters);

protected:
  virtual ADReal computeQpResidual() override;

  // A vector to store all of the isotope number densities which form the
  // current isotope under neutron bombardment.
  std::vector<const ADVariableValue *> _isotope_fractions;

  // A vector to store all of the decay constants.
  std::vector<Real> _decay_consts;
}; // class ADMassFractionNuclideDecaySource
