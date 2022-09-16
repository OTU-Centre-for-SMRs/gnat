#pragma once

#include "ADIsotopeBaseBC.h"

class ADIsotopeInflowBC : public ADIsotopeBaseBC
{
public:
  static InputParameters validParams();

  ADIsotopeInflowBC(const InputParameters & parameters);

protected:
  virtual ADReal computeQpResidual() override;

  const Real _inflow_rate;
};
