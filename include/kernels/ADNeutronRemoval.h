#pragma once

#include "ADKernel.h"

class ADNeutronRemoval : public ADKernel
{
public:
  static InputParameters validParams();

  ADNeutronRemoval(const InputParameters & parameters);

protected:
  virtual ADReal computeQpResidual() override;

  const unsigned int _group_index; // g

  const ADMaterialProperty<std::vector<Real>> & _sigma_r_g;
};
