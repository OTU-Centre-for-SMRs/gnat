#pragma once

#include "ADIsotopeBase.h"

class ADIsotopeAdvection : public ADIsotopeBase
{
public:
  static InputParameters validParams();

  ADIsotopeAdvection(const InputParameters & parameters);

protected:
  virtual ADReal computeQpResidual() override;
}; // class ADIsotopeAdvection
