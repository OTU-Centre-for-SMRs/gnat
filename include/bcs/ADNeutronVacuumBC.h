#pragma once

#include "ADIntegratedBC.h"

class ADNeutronVacuumBC : public ADIntegratedBC
{
public:
  static InputParameters validParams();

  ADNeutronVacuumBC(const InputParameters & parameters);

protected:
  virtual ADReal computeQpResidual() override;

  const unsigned int _ordinate_index; // n

  const MaterialProperty<std::vector<RealVectorValue>> & _directions;
}; // class ADNeutronVacuumBC
