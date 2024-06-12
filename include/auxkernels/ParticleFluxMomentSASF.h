#pragma once

#include "AuxKernel.h"

// An auxkernel which computes the other angular flux moments from a point source given the scalar
// flux.
class ParticleFluxMomentSASF : public AuxKernel
{
public:
  static InputParameters validParams();

  ParticleFluxMomentSASF(const InputParameters & parameters);

protected:
  virtual Real computeValue() override;
  void cartesianToSpherical(const RealVectorValue & direction, Real & mu, Real & omega);

  // The source location.
  const Point _source_location;

  const unsigned int _degree;
  const int _order;

  const VariableValue & _uncollided_scalar_flux;
}; // class ParticleFluxMomentSASF
