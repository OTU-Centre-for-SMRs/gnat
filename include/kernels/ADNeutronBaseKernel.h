#pragma once

#include "ADKernel.h"

#include "GaussAngularQuadrature.h"

// A base class for all neutron transport kernels which provides basic
// components that appear in the terms of the SN neutron transport equation weak
// form.
class ADNeutronBaseKernel : public ADKernel
{
public:
  static InputParameters validParams();

  ADNeutronBaseKernel(const InputParameters & parameters);

protected:
  void cartesianToSpherical(const RealVectorValue & ordinate,
                            Real & mu, Real & omega);

  const GaussAngularQuadrature _quadrature_set;
  Real _symmetry_factor;
}; // class ADNeutronBaseKernel
