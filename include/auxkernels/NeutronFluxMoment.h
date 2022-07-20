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

  void cartesianToSpherical(const RealVectorValue & ordinate, Real & mu,
                            Real & omega);

  const GaussAngularQuadrature _quadrature_set;

  std::vector<const ADVariableValue *> _flux_ordinates;

  const unsigned int _degree;
  const int _order;
}; // class NeutronFluxMoment
