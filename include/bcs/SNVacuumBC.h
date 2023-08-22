#pragma once

#include "SNBaseBC.h"

class SNVacuumBC : public SNBaseBC
{
public:
  static InputParameters validParams();

  SNVacuumBC(const InputParameters & parameters);

protected:
  virtual Real computeQpResidual() override;
  virtual Real computeQpJacobian() override;

  const unsigned int _ordinate_index; // n
};                                    // class SNVacuumBC
