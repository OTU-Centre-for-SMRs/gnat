#pragma once

#include "ADIntegratedBC.h"

#include "GaussAngularQuadrature.h"

class ADNeutronBaseBC : public ADIntegratedBC
{
public:
  static InputParameters validParams();

  ADNeutronBaseBC(const InputParameters & parameters);

protected:
  void cartesianToSpherical(const RealVectorValue & ordinate,
                            Real & mu, Real & omega);

  const GaussAngularQuadrature _quadrature_set;
  Real _symmetry_factor;
}; // class ADNeutronBaseBC
