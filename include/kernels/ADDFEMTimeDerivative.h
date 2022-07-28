#pragma once

#include "ADTimeDerivative.h"

class ADDFEMTimeDerivative : public ADTimeDerivative
{
public:
  static InputParameters validParams();

  ADDFEMTimeDerivative(const InputParameters & parameters);
protected:
  virtual ADReal precomputeQpResidual() override;

  const unsigned int _group_index; // g

  const ADMaterialProperty<std::vector<Real>> & _v_g;
}; // class ADDFEMTimeDerivative
