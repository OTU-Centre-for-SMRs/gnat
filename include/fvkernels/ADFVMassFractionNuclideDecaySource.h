#pragma once

#include "FVElementalKernel.h"

class ADFVMassFractionNuclideDecaySource : public FVElementalKernel
{
public:
  static InputParameters validParams();

  ADFVMassFractionNuclideDecaySource(const InputParameters & parameters);

protected:
  virtual ADReal computeQpResidual() override;

  // Density of the bulk fluid.
  const Moose::Functor<ADReal> & _density;

  // A vector to store all of the isotope number densities which form the
  // current isotope when decaying..
  std::vector<const Moose::Functor<ADReal> *> _isotope_fractions;

  // A vector to store all of the decay constants.
  std::vector<Real> _decay_consts;
}; // class ADFVMassFractionNuclideDecaySource
