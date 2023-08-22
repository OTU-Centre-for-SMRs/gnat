#pragma once

#include "ADIsotopeBaseBC.h"

class ADMassFractionNuclideOutflowBC : public ADIsotopeBaseBC
{
public:
  static InputParameters validParams();

  ADMassFractionNuclideOutflowBC(const InputParameters & parameters);

protected:
  virtual ADReal computeQpResidual() override;
}; // class ADMassFractionNuclideOutflowBC
