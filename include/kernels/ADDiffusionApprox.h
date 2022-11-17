#pragma once

#include "ADKernel.h"

// Kernel to compute the residual contribution for the diffusion approximation in the neutron
// transport equation.
class ADDiffusionApprox : public ADKernel
{
public:
  static InputParameters validParams();

  ADDiffusionApprox(const InputParameters & parameters);

protected:
  virtual ADReal computeQpResidual() override;

  // g
  const unsigned int _group_index;

  const ADMaterialProperty<std::vector<Real>> & _diffusion_g;
}; // class ADDiffusionApprox
