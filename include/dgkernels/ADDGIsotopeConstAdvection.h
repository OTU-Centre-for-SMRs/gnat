#pragma once

#include "ADDGKernel.h"

class ADDGIsotopeConstAdvection : public ADDGKernel
{
public:
  static InputParameters validParams();

  ADDGIsotopeConstAdvection(const InputParameters & parameters);

protected:
  virtual ADReal computeQpResidual(Moose::DGResidualType type) override;

  const RealVectorValue _vel;
};
