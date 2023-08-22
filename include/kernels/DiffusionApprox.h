#pragma once

#include "Kernel.h"

// Kernel to compute the residual contribution for the diffusion approximation in the neutron
// transport equation.
class DiffusionApprox : public Kernel
{
public:
  static InputParameters validParams();

  DiffusionApprox(const InputParameters & parameters);

protected:
  virtual Real computeQpResidual() override;
  virtual Real computeQpJacobian() override;

  // g
  const unsigned int _group_index;

  const ADMaterialProperty<std::vector<Real>> & _diffusion_g;
}; // class DiffusionApprox
