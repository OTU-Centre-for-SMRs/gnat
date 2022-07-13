#pragma once

#include "ADKernel.h"

class ADNeutronStreaming : public ADKernel
{
public:
  static InputParameters validParams();

  ADNeutronStreaming(const InputParameters & parameters);

protected:
  virtual ADReal computeQpResidual() override;

  const unsigned int _ordinate_index; // n

  const ADMaterialProperty<std::vector<RealVectorValue>> & _directions;
}; // class ADNeutronStreaming
