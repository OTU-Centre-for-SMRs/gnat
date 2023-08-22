#pragma once

#include "ADIsotopeBaseBC.h"

class ADMassFractionNuclideInflowBC : public ADIsotopeBaseBC
{
public:
  static InputParameters validParams();

  ADMassFractionNuclideInflowBC(const InputParameters & parameters);

protected:
  virtual ADReal computeQpResidual() override;

  const Real _inflow_rate;
}; // class ADMassFractionNuclideInflowBC
