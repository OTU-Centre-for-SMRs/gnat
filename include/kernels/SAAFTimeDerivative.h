#pragma once

#include "SAAFBaseKernel.h"

class SAAFTimeDerivative : public SAAFBaseKernel
{
public:
  static InputParameters validParams();

  SAAFTimeDerivative(const InputParameters & parameters);

protected:
  virtual Real computeQpResidual() override;
  virtual Real computeQpJacobian() override;

  // Holds the time derivatives at the quadrature points.
  const VariableValue & _u_dot;
  const VariableValue & _du_dot_du;

  const ADMaterialProperty<std::vector<Real>> & _inv_v_g;
}; // ADSAAFTimeDerivative
