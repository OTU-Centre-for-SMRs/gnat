#pragma once

#include "ADIntegratedBC.h"

class ADNeutronMatSourceBC : public ADIntegratedBC
{
public:
  static InputParameters validParams();

  ADNeutronMatSourceBC(const InputParameters & parameters);

protected:
  virtual ADReal computeQpResidual() override;

  const unsigned int _ordinate_index; // n

  const ADMaterialProperty<std::vector<RealVectorValue>> & _directions;
  const ADMaterialProperty<std::vector<Real>> & _surface_source;
}; // class ADNeutronMatSourceBC
