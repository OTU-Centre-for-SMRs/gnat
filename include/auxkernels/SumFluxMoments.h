#pragma once

#include "AuxKernel.h"

// An auxkernel to compute the sum of the collided and uncollided components of flux moments.
class SumFluxMoments : public AuxKernel
{
public:
  static InputParameters validParams();

  SumFluxMoments(const InputParameters & parameters);

protected:
  virtual Real computeValue() override;

  // Collided and uncollided moments of the angular particle flux.
  const ADVariableValue & _collided_moment;
  const ADVariableValue & _uncollided_moment;
}; // class SumFluxMoments
