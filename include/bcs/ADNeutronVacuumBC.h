#pragma once

#include "ADNeutronBaseBC.h"

class ADNeutronVacuumBC : public ADNeutronBaseBC
{
public:
  static InputParameters validParams();

  ADNeutronVacuumBC(const InputParameters & parameters);

protected:
  virtual ADReal computeQpResidual() override;

  const unsigned int _ordinate_index; // n
}; // class ADNeutronVacuumBC
