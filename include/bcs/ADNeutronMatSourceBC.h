#pragma once

#include "ADNeutronBaseBC.h"

class ADNeutronMatSourceBC : public ADNeutronBaseBC
{
public:
  static InputParameters validParams();

  ADNeutronMatSourceBC(const InputParameters & parameters);

protected:
  virtual ADReal computeQpResidual() override;

  const unsigned int _ordinate_index; // n

  const ADMaterialProperty<std::vector<Real>> & _surface_source;
}; // class ADNeutronMatSourceBC
