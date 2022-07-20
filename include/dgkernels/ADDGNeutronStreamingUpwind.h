#pragma once

#include "ADDGKernel.h"

#include "GaussAngularQuadrature.h"

class ADDGNeutronStreamingUpwind : public ADDGKernel
{
public:
  static InputParameters validParams();

  ADDGNeutronStreamingUpwind(const InputParameters & parameters);

protected:
  virtual ADReal computeQpResidual(Moose::DGResidualType type) override;

  const GaussAngularQuadrature _quadrature_set;

  const unsigned int _ordinate_index; // n
}; // class ADDGNeutronStreamingUpwind
