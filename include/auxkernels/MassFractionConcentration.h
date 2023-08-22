#pragma once

#include "AuxKernel.h"

class MassFractionConcentration : public AuxKernel
{
public:
  static InputParameters validParams();

  MassFractionConcentration(const InputParameters & parameters);

protected:
  static constexpr Real _n_avogadro = 6.0221408e23;

  virtual Real computeValue() override;

  const bool _output_number_density;

  // Density of the bulk fluid.
  const Moose::Functor<ADReal> & _density;

  // The mass fraction variable in either variable (FE) or functor (FV) form.
  const bool _is_fe;
  const ADVariableValue * _mass_fraction_var;
  const Moose::Functor<ADReal> * _mass_fraction_fun;

  // The atomic weight (molar mass) of the nuclide.
  Real _weight;
};
