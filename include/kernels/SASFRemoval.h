#pragma once

#include "Kernel.h"

class SASFRemoval : public Kernel
{
public:
  static InputParameters validParams();

  SASFRemoval(const InputParameters & parameters);

protected:
  virtual Real computeQpResidual() override;
  virtual Real computeQpJacobian() override;

  // g
  const unsigned int _group_index;
  // Total cross-section.
  const ADMaterialProperty<std::vector<Real>> & _sigma_t_g;
}; // class SASFRemoval
