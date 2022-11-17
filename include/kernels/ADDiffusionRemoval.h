#pragma once

#include "ADKernel.h"

class ADDiffusionRemoval : public ADKernel
{
public:
  static InputParameters validParams();

  ADDiffusionRemoval(const InputParameters & parameters);

protected:
  virtual ADReal computeQpResidual() override;

  // g
  const unsigned int _group_index;

  const ADMaterialProperty<std::vector<Real>> & _sigma_r_g;
}; // class ADDiffusionRemoval
