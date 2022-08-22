#pragma once

#include "ADKernel.h"

class ADIsotopeDiffusion : public ADKernel
{
public:
  static InputParameters validParams();

  ADIsotopeDiffusion(const InputParameters & parameters);

protected:
  virtual ADReal computeQpResidual() override;

  const ADMaterialProperty<Real> & _mat_diff;
};
