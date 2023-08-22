#pragma once

#include "ADIsotopeBase.h"

class ADMassFractionNuclideAdvection : public ADIsotopeBase
{
public:
  static InputParameters validParams();

  ADMassFractionNuclideAdvection(const InputParameters & parameters);

protected:
  virtual ADReal computeQpResidual() override;

  ADReal computeQpVelDivergence();
}; // class ADMassFractionNuclideAdvection
