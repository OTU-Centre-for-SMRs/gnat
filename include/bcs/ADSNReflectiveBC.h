#pragma once

#include "ADSNBaseBC.h"

class ADSNReflectiveBC : public ADSNBaseBC
{
public:
  static InputParameters validParams();

  ADSNReflectiveBC(const InputParameters & parameters);

protected:
  virtual ADReal computeQpResidual() override;

  const unsigned int _ordinate_index; // n

  const ADVariableValue & _u_ref;

  const ADMaterialProperty<Real> & _albedo;
}; // class ADSNReflectiveBC
