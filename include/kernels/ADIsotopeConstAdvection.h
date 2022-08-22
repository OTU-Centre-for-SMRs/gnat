#pragma once

#include "ADKernel.h"

class ADIsotopeConstAdvection : public ADKernel
{
public:
  static InputParameters validParams();

  ADIsotopeConstAdvection(const InputParameters & parameters);

protected:
  virtual ADReal computeQpResidual() override;

  const RealVectorValue _vel;
};
