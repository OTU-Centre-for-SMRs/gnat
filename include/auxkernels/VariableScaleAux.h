#pragma once

#include "AuxKernel.h"

class VariableScaleAux : public AuxKernel
{
public:
  static InputParameters validParams();

  VariableScaleAux(const InputParameters & parameters);

protected:
  virtual Real computeValue() override;

  // The concentration variable in either variable (FE) or functor (FV) form.
  const bool _is_fe;
  const ADVariableValue * _to_scale_var;
  const Moose::Functor<ADReal> * _to_scale_fun;

  // The value to multiply the field by.
  const Real _scale_factor;
}; // class VariableScaleAux
