#pragma once

#include "ADKernel.h"

class ADSNRemoval : public ADKernel
{
public:
  static InputParameters validParams();

  ADSNRemoval(const InputParameters & parameters);

protected:
  virtual ADReal computeQpResidual() override;

  // g
  const unsigned int _group_index;

  const ADMaterialProperty<std::vector<Real>> & _sigma_r_g;
}; // class ADSNRemoval
