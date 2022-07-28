#pragma once

#include "ADIntegratedBC.h"

#include "GaussAngularQuadrature.h"

class ADSNBaseBC : public ADIntegratedBC
{
public:
  static InputParameters validParams();

  ADSNBaseBC(const InputParameters & parameters);

protected:
  void cartesianToSpherical(const RealVectorValue & ordinate,
                            Real & mu, Real & omega);

  const GaussAngularQuadrature _quadrature_set;
  Real _symmetry_factor;
}; // class ADSNBaseBC
