#pragma once

#include "ADSNBaseBC.h"

class ADSNMatSourceBC : public ADSNBaseBC
{
public:
  static InputParameters validParams();

  ADSNMatSourceBC(const InputParameters & parameters);

protected:
  virtual ADReal computeQpResidual() override;

  const unsigned int _ordinate_index; // n

  const ADMaterialProperty<std::vector<Real>> & _surface_source;
}; // class ADSNMatSourceBC
