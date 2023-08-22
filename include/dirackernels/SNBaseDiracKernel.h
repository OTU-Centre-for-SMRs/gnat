#pragma once

#include "DiracKernel.h"

#include "AQProvider.h"

class SNBaseDiracKernel : public DiracKernel
{
public:
  static InputParameters validParams();

  SNBaseDiracKernel(const InputParameters & parameters);

protected:
  void cartesianToSpherical(const RealVectorValue & ordinate, Real & mu,
                            Real & omega);

  const AQProvider & _aq;
  Real _symmetry_factor;
}; // class SNBaseDiracKernel
