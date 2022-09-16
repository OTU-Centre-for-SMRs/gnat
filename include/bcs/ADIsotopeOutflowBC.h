#pragma once

#include "ADIsotopeBaseBC.h"

class ADIsotopeOutflowBC : public ADIsotopeBaseBC
{
public:
  static InputParameters validParams();

  ADIsotopeOutflowBC(const InputParameters & parameters);

protected:
  virtual ADReal computeQpResidual() override;
};
