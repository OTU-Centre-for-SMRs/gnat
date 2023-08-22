#pragma once

#include "IntegratedBC.h"

#include "AQProvider.h"

class SNBaseBC : public IntegratedBC
{
public:
  static InputParameters validParams();

  SNBaseBC(const InputParameters & parameters);

protected:
  void cartesianToSpherical(const RealVectorValue & ordinate, Real & mu, Real & omega);

  const AQProvider & _aq;
  Real _symmetry_factor;
}; // class SNBaseBC
