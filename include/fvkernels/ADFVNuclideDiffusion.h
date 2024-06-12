#pragma once

#include "FVFluxKernel.h"

class ADFVNuclideDiffusion : public FVFluxKernel
{
public:
  static InputParameters validParams();

  ADFVNuclideDiffusion(const InputParameters & params);

protected:
  ADReal computeQpResidual() override;

  const Moose::Functor<ADReal> & _mat_diff;
}; // class ADFVNuclideDiffusion
