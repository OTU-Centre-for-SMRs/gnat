#pragma once

#include "ADKernel.h"

// TODO: Finish implementing this kernel.
class ADNeutronGToGScattering : public ADKernel
{
public:
  static InputParameters validParams();

  ADNeutronGToGScattering(const InputParameters & parameters);

protected:
  virtual ADReal computeQpResidual() override;
};
