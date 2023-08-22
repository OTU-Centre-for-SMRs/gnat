#pragma once

#include "ADIsotopeBase.h"

class ADMassFractionNuclideTimeDerivative : public ADIsotopeBase
{
public:
  static InputParameters validParams();

  ADMassFractionNuclideTimeDerivative(const InputParameters & parameters);

protected:
  virtual ADReal computeQpResidual() override;

  // Holds the time derivatives at the quadrature points.
  const ADTemplateVariableValue<Real> & _u_dot;
}; // class ADMassFractionNuclideTimeDerivative
