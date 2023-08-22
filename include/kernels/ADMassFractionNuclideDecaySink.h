#pragma once

#include "ADIsotopeBase.h"

class ADMassFractionNuclideDecaySink : public ADIsotopeBase
{
public:
  static InputParameters validParams();

  ADMassFractionNuclideDecaySink(const InputParameters & parameters);

protected:
  virtual ADReal computeQpResidual() override;

  const Real _decay_const;
}; // class ADMassFractionNuclideDecaySink
