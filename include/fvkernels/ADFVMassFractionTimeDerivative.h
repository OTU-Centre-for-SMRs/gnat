#pragma once

#include "FVTimeKernel.h"

class ADFVMassFractionTimeDerivative : public FVTimeKernel
{
public:
  static InputParameters validParams();

  ADFVMassFractionTimeDerivative(const InputParameters & parameters);

protected:
  ADReal computeQpResidual() override;

  // Density of the bulk fluid.
  const Moose::Functor<ADReal> & _density;
};
