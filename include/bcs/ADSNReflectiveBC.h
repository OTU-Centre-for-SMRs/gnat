#pragma once

#include "ADSNBaseBC.h"

class ADSNReflectiveBC : public ADSNBaseBC
{
public:
  static InputParameters validParams();

  ADSNReflectiveBC(const InputParameters & parameters);

protected:
  virtual ADReal computeQpResidual() override;
  unsigned int findReflectedOrdinate();

  const unsigned int _ordinate_index; // n

  // All angular fluxes in the current group.
  std::vector<const VariableValue *> _flux_ordinates;
}; // class ADSNReflectiveBC
