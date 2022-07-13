#pragma once

#include "ADIntegratedBC.h"

class ADNeutronReflectiveBC : public ADIntegratedBC
{
public:
  static InputParameters validParams();

  ADNeutronReflectiveBC(const InputParameters & parameters);
protected:
  virtual ADReal computeQpResidual() override;

  const unsigned int _ordinate_index; // n

  const ADVariableValue & _u_ref;

  const ADMaterialProperty<std::vector<RealVectorValue>> & _directions;
  const ADMaterialProperty<Real> & _albedo;
}; // class ADNeutronReflectiveBC
