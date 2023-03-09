#pragma once

#include "SNBaseKernel.h"

class SAAFBaseKernel : public SNBaseKernel
{
public:
  static InputParameters validParams();

  SAAFBaseKernel(const InputParameters & parameters);

protected:
  // Computes $\tau$ for the current quadrature point.
  Real computeQPTau();
  // Computes $\phi_{j} + \tau_{g}\vec{\nabla}\phi_{j}\cdot\hat{\Omega}$ for the
  // current quadrature point.
  Real computeQPTests();

  // n
  const unsigned int _ordinate_index;
  // g
  const unsigned int _group_index;

  const ADMaterialProperty<std::vector<Real>> & _sigma_r_g;

  // SAAF stabilization parameters.
  const ADMaterialProperty<Real> & _saaf_eta;
  const ADMaterialProperty<Real> & _saaf_c;
}; // class SAAFBaseKernel
