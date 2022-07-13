#pragma once

#include "ADDGKernel.h"

class ADDGNeutronStreamingUpwind : public ADDGKernel
{
public:
  static InputParameters validParams();

  ADDGNeutronStreamingUpwind(const InputParameters & parameters);

protected:
  virtual ADReal computeQpResidual(Moose::DGResidualType type) override;

  const unsigned int _ordinate_index; // n

  const ADMaterialProperty<std::vector<RealVectorValue>> & _directions;
}; // class ADDGNeutronStreamingUpwind
