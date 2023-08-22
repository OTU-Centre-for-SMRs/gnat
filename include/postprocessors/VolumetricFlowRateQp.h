#pragma once

#include "SideIntegralPostprocessor.h"

// A modified version of the 'VolumetricFlowRate' postprocessor from the MOOSE Navier-Stokes module.
// This couples finite volume velocity fields to finite element advected quantities through
// quadrature points.
class VolumetricFlowRateQp : public SideIntegralPostprocessor
{
public:
  static InputParameters validParams();

  VolumetricFlowRateQp(const InputParameters & parameters);

protected:
  Real computeQpIntegral() override;

  const unsigned int _mesh_dims;

  // The velocity field.
  // X-velocity.
  const Moose::Functor<ADReal> & _vel_u;
  // Y-velocity.
  const Moose::Functor<ADReal> * const _vel_v;
  // Z-velocity.
  const Moose::Functor<ADReal> * const _vel_w;

  // The advected quantity.
  const Moose::Functor<ADReal> & _adv_quant;
}; // class VolumetricFlowRateQp
