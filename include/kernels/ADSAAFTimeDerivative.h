#pragma once

#include "ADSAAFBaseKernel.h"

class ADSAAFTimeDerivative : public ADSAAFBaseKernel
{
public:
  static InputParameters validParams();

  ADSAAFTimeDerivative(const InputParameters & parameters);

protected:
  virtual ADReal computeQpResidual() override;

  // Holds the time derivatives at the quadrature points.
  const ADTemplateVariableValue<Real> & _u_dot;

  const ADMaterialProperty<std::vector<Real>> & _inv_v_g;
}; // ADSAAFTimeDerivative
