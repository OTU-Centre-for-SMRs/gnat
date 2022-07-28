#pragma once

#include "ADDGKernel.h"

#include "GaussAngularQuadrature.h"

class ADDFEMUpwinding : public ADDGKernel
{
public:
  static InputParameters validParams();

  ADDFEMUpwinding(const InputParameters & parameters);

protected:
  virtual ADReal computeQpResidual(Moose::DGResidualType type) override;

  const GaussAngularQuadrature _quadrature_set;

  const unsigned int _ordinate_index; // n
}; // class ADDFEMUpwinding
