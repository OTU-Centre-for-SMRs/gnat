#pragma once

#include "ADIsotopeBase.h"

class ADIsotopeDiffusion : public ADIsotopeBase
{
public:
  static InputParameters validParams();

  ADIsotopeDiffusion(const InputParameters & parameters);

protected:
  virtual ADReal computeQpResidual() override;

  const ADMaterialProperty<Real> & _mat_diff;
  const ADMaterialProperty<RealVectorValue> & _grad_mat_diff;
}; // class ADIsotopeDiffusion
