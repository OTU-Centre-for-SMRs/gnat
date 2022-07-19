#pragma once

#include "ADIntegratedBC.h"
#include "GnatBase.h"

class ADNeutronMatSourceBC : public ADIntegratedBC
{
public:
  static InputParameters validParams();

  ADNeutronMatSourceBC(const InputParameters & parameters);

protected:
  virtual ADReal computeQpResidual() override;

  const ProblemType _type;
  Real _symmetry_factor;

  const unsigned int _ordinate_index; // n

  const MaterialProperty<std::vector<RealVectorValue>> & _directions;
  const ADMaterialProperty<std::vector<Real>> & _surface_source;
}; // class ADNeutronMatSourceBC
