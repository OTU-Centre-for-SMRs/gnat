#pragma once

#include "FVElementalKernel.h"

class ADFVNuclideDecaySource : public FVElementalKernel
{
public:
  static InputParameters validParams();

  ADFVNuclideDecaySource(const InputParameters & parameters);

protected:
  virtual ADReal computeQpResidual() override;

  // A vector to store all of the isotope number densities which form the
  // current isotope when decaying..
  std::vector<const Moose::Functor<ADReal> *> _isotope_number_densities;

  // A vector to store all of the decay constants.
  std::vector<Real> _decay_consts;
}; // class ADFVNuclideDecaySource
