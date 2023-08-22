#pragma once

#include "FVFluxKernel.h"

class ADFVMassFractionNuclideDiffusion : public FVFluxKernel
{
public:
  static InputParameters validParams();

  ADFVMassFractionNuclideDiffusion(const InputParameters & params);

protected:
  ADReal computeQpResidual() override;

  // Density of the bulk fluid.
  const Moose::Functor<ADReal> & _density;

  const Moose::Functor<ADReal> & _mat_diff;
}; // class ADFVMassFractionNuclideDiffusion
