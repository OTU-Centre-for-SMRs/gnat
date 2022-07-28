#pragma once

#include "DiracKernel.h"

#include "GaussAngularQuadrature.h"

class SNDiracKernelBase : public DiracKernel
{
public:
  static InputParameters validParams();

  SNDiracKernelBase(const InputParameters & parameters);

protected:
  void cartesianToSpherical(const RealVectorValue & ordinate, Real & mu,
                            Real & omega);

  const GaussAngularQuadrature _quadrature_set;
  Real _symmetry_factor;
}; // class SNDiracKernelBase
