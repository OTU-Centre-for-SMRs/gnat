#pragma once

#include "ADTimeDerivative.h"

class ADDFEMTimeDerivative : public ADTimeDerivative
{
public:
  static InputParameters validParams();

  ADDFEMTimeDerivative(const InputParameters & parameters);

protected:
  virtual ADReal precomputeQpResidual() override;

  // g
  const unsigned int _group_index;

  const ADMaterialProperty<std::vector<Real>> & _v_g;
}; // class ADDFEMTimeDerivative
