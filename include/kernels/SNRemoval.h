#pragma once

#include "Kernel.h"

class SNRemoval : public Kernel
{
public:
  static InputParameters validParams();

  SNRemoval(const InputParameters & parameters);

protected:
  virtual Real computeQpResidual() override;
  virtual Real computeQpJacobian() override;

  // g
  const unsigned int _group_index;

  const ADMaterialProperty<std::vector<Real>> & _sigma_t_g;
}; // class ADSNRemoval
