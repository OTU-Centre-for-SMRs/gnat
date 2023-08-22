#pragma once

#include "AuxKernel.h"

#include "AQProvider.h"

class ParticleFluxMoment : public AuxKernel
{
public:
  static InputParameters validParams();

  ParticleFluxMoment(const InputParameters & parameters);

protected:
  virtual Real computeValue() override;

  const AQProvider & _aq;

  std::vector<const ADVariableValue *> _flux_ordinates;

  const unsigned int _degree;
  const int _order;

  const Real _scale_factor;
}; // class ParticleFluxMoment
