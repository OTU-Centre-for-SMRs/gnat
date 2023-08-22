#pragma once

#include "AuxKernel.h"

// A class to compute the specific activity (Bq / cm^{-3}) of a radionuclide.
class SpecificActivity : public AuxKernel
{
public:
  static InputParameters validParams();

  SpecificActivity(const InputParameters & parameters);

protected:
  virtual Real computeValue() override;

  // The concentration variable in either variable (FE) or functor (FV) form.
  const bool _is_fe;
  const ADVariableValue * _concentration_var;
  const Moose::Functor<ADReal> * _concentration_fun;

  // Decay constant.
  Real _decay_const;
}; // class SpecificActivity
