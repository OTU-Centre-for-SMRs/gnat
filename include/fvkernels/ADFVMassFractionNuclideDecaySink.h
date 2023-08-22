#pragma once

#include "FVElementalKernel.h"

class ADFVMassFractionNuclideDecaySink : public FVElementalKernel
{
public:
  static InputParameters validParams();

  ADFVMassFractionNuclideDecaySink(const InputParameters & parameters);

protected:
  virtual ADReal computeQpResidual() override;

  // Density of the bulk fluid.
  const Moose::Functor<ADReal> & _density;

  const Real _decay_const;
}; // class ADFVMassFractionNuclideDecaySink
