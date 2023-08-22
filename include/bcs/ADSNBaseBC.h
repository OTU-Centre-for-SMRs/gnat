#pragma once

#include "ADIntegratedBC.h"

#include "AQProvider.h"

class ADSNBaseBC : public ADIntegratedBC
{
public:
  static InputParameters validParams();

  ADSNBaseBC(const InputParameters & parameters);

protected:
  void cartesianToSpherical(const RealVectorValue & ordinate, Real & mu, Real & omega);

  const AQProvider & _aq;
  Real _symmetry_factor;
}; // class ADSNBaseBC
