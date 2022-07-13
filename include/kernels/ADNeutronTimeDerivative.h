#pragma once

#include "ADTimeDerivative.h"

class ADNeutronTimeDerivative : public ADTimeDerivative
{
public:
  static InputParameters validParams();

  ADNeutronTimeDerivative(const InputParameters & parameters);
protected:
  virtual ADReal precomputeQpResidual() override;

  const unsigned int _group_index; // g

  const ADMaterialProperty<std::vector<Real>> & _v_g;
}; // class ADNeutronTimeDerivative
