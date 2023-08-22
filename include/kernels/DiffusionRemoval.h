#pragma once

#include "Kernel.h"

class DiffusionRemoval : public Kernel
{
public:
  static InputParameters validParams();

  DiffusionRemoval(const InputParameters & parameters);

protected:
  virtual Real computeQpResidual() override;
  virtual Real computeQpJacobian() override;

  // g
  const unsigned int _group_index;

  const ADMaterialProperty<std::vector<Real>> & _sigma_r_g;
}; // class DiffusionRemoval
