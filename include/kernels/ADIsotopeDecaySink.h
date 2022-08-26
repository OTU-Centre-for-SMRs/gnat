#pragma once

#include "ADIsotopeBase.h"

class ADIsotopeDecaySink : public ADIsotopeBase
{
public:
  static InputParameters validParams();

  ADIsotopeDecaySink(const InputParameters & parameters);

protected:
  virtual ADReal computeQpResidual() override;

  const Real _decay_const;
}; // class ADIsotopeDecaySink
