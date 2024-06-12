#pragma once

#include "FVElementalKernel.h"

class ADFVNuclideDecaySink : public FVElementalKernel
{
public:
  static InputParameters validParams();

  ADFVNuclideDecaySink(const InputParameters & parameters);

protected:
  virtual ADReal computeQpResidual() override;

  const Real _decay_const;
}; // class ADFVNuclideDecaySink
