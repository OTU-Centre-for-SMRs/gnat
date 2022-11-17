#pragma once

#include "ADTimeDerivative.h"

class ADParticleTimeDerivative : public ADTimeDerivative
{
public:
  static InputParameters validParams();

  ADParticleTimeDerivative(const InputParameters & parameters);

protected:
  virtual ADReal precomputeQpResidual() override;

  // g
  const unsigned int _group_index;

  const ADMaterialProperty<std::vector<Real>> & _inv_v_g;
}; // class ADParticleTimeDerivative
