#pragma once

#include "ADKernel.h"

class ADIsotopeDecaySink : public ADKernel
{
public:
  static InputParameters validParams();

  ADIsotopeDecaySink(const InputParameters & parameters);

protected:
  virtual ADReal computeQpResidual() override;

  const Real _decay_const;
}; // class ADIsotopeDecaySink
