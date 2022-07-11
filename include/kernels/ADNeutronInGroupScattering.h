#pragma once

#include "ADKernel.h"

// TODO: Finish implementing this kernel.
class ADNeutronInGroupScattering : public ADKernel
{
public:
  static InputParameters validParams();

  ADNeutronInGroupScattering(const InputParameters & parameters);

protected:
  virtual ADReal computeQpResidual() override;
};
