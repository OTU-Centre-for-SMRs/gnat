#pragma once

#include "INSFVAdvectionKernel.h"

// A modified version of INSFVScalarFieldAdvection to handle advecting a trace species which is
// defined based on a mass fraction of the bulk fluid.
class INSFVMassFractionScalarFieldAdvection : public INSFVAdvectionKernel
{
public:
  static InputParameters validParams();
  INSFVMassFractionScalarFieldAdvection(const InputParameters & params);

protected:
  ADReal computeQpResidual() override;

  virtual bool hasMaterialTimeDerivative() const override { return true; }

  // Density of the bulk fluid.
  const Moose::Functor<ADReal> & _density;
}; // class INSFVMassFractionScalarFieldAdvection
