#pragma once

#include "ADKernel.h"

class ADSNRemoval : public ADKernel
{
public:
  static InputParameters validParams();

  ADSNRemoval(const InputParameters & parameters);

protected:
  virtual ADReal computeQpResidual() override;

  const unsigned int _group_index; // g

  const ADMaterialProperty<std::vector<Real>> & _sigma_r_g;
}; // class ADSNRemoval
