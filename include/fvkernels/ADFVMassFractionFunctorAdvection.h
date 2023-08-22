#pragma once

#include "FVFluxKernel.h"

#include "MathFVUtils.h"

class ADFVMassFractionFunctorAdvection : public FVFluxKernel
{
public:
  static InputParameters validParams();

  ADFVMassFractionFunctorAdvection(const InputParameters & parameters);

protected:
  ADReal computeQpResidual() override;

  const unsigned int _mesh_dims;

  // Density of the bulk fluid.
  const Moose::Functor<ADReal> & _density;

  const Moose::Functor<ADReal> & _vel_u;
  // Y-velocity.
  const Moose::Functor<ADReal> * const _vel_v;
  // Z-velocity.
  const Moose::Functor<ADReal> * const _vel_w;

  /// The interpolation method to use for the advected quantity
  Moose::FV::InterpMethod _advected_interp_method;

  /// The interpolation method to use for the velocity
  Moose::FV::InterpMethod _velocity_interp_method;
}; // class ADFVMassFractionFunctorAdvection
