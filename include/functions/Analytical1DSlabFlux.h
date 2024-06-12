#pragma once

#include "Function.h"

// A function which computes the analytical scalar flux for an isotropic planar source embedded in a
// 1D infinite absorbing slab
class Analytical1DSlabFlux : public Function
{
public:
  static InputParameters validParams();

  Analytical1DSlabFlux(const InputParameters & parameters);

  virtual Real value(Real t, const Point & p) const override;

protected:
  // The planar source emission strength.
  const Real _source_strength;
  // The location of the planar source.
  const Real _source_location;

  // The total macroscopic cross-section of the slab.
  const Real & _total_xs;
}; // class Analytical1DSlabFlux.
