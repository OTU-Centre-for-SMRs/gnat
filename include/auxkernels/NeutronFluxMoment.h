#pragma once

#include "AuxKernel.h"

#include "GaussAngularQuadrature.h"

class NeutronFluxMoment : public AuxKernel
{
public:
  static InputParameters validParams();

  NeutronFluxMoment(const InputParameters & parameters);

protected:
  virtual Real computeValue() override;

  void cartesianToSpherical(const RealVectorValue & ordinate, Real & mu, Real & omega);

  const GaussAngularQuadrature _quadrature_set;

  const enum class NegativeFluxHandling { Max = 0u, Abs = 1u, None = 2u } _flux_handling;

  std::vector<const ADVariableValue *> _flux_ordinates;

  const unsigned int _degree;
  const int _order;
  Real _symmetry_factor;
}; // class NeutronFluxMoment
