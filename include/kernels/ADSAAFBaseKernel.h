#pragma once

#include "ADSNBaseKernel.h"

class ADSAAFBaseKernel : public ADSNBaseKernel
{
public:
  static InputParameters validParams();

  ADSAAFBaseKernel(const InputParameters & parameters);

protected:
  // Computes $\tau$ for the current quadrature point.
  ADReal computeQPTau();
  // Computes $\phi_{j} + \tau_{g}\vec{\nabla}\phi_{j}\cdot\hat{\Omega}$ for the
  // current quadrature point.
  ADReal computeQPTests();

  const unsigned int _ordinate_index; // n
  const unsigned int _group_index; // g

  const ADMaterialProperty<std::vector<Real>> & _sigma_r_g;

  // SAAF stabilization parameters.
  const ADMaterialProperty<Real> & _saaf_eta;
  const ADMaterialProperty<Real> & _saaf_c;
}; // class ADSAAFBaseKernel
