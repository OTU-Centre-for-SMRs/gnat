#pragma once

#include "IntegratedBC.h"

class SASFVacuumBC : public IntegratedBC
{
public:
  static InputParameters validParams();

  SASFVacuumBC(const InputParameters & parameters);

protected:
  virtual Real computeQpResidual() override;
  virtual Real computeQpJacobian() override;

  // The source location.
  const Point _source_location;
}; // class SASFVacuumBC
