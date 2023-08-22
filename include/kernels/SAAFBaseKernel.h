#pragma once

#include "SNBaseKernel.h"

// A class which provides common functionality to discrete ordinates kernels which have been
// stabilized with the self-adjoint angular flux method.
class SAAFBaseKernel : public SNBaseKernel
{
public:
  static InputParameters validParams();

  SAAFBaseKernel(const InputParameters & parameters);

protected:
  // Computes $\phi_{j} + \tau_{g}\vec{\nabla}\phi_{j}\cdot\hat{\Omega}$ for the
  // current quadrature point.
  Real computeQpTests();

  // n
  const unsigned int _ordinate_index;
  // g
  const unsigned int _group_index;

  const ADMaterialProperty<std::vector<Real>> & _sigma_t_g;

  // SAAF stabilization parameters.
  const ADMaterialProperty<std::vector<Real>> & _saaf_tau;
}; // class SAAFBaseKernel
