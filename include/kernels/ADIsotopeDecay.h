#pragma once

#include "ADKernel.h"

class ADIsotopeDecay : public ADKernel
{
public:
  static InputParameters validParams();

  ADIsotopeDecay(const InputParameters & parameters);

protected:
  virtual ADReal computeQpResidual() override;

  const Real _decay_const;
};
