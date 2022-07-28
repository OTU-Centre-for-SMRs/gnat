#pragma once

#include "ADSNBaseBC.h"

class ADSNVacuumBC : public ADSNBaseBC
{
public:
  static InputParameters validParams();

  ADSNVacuumBC(const InputParameters & parameters);

protected:
  virtual ADReal computeQpResidual() override;

  const unsigned int _ordinate_index; // n
}; // class ADSNVacuumBC
