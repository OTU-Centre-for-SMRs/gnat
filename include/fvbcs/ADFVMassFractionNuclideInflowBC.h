#pragma once

#include "FVFluxBC.h"

// A class to impose a concentration along a boundary given a flux.
class ADFVMassFractionNuclideInflowBC : public FVFluxBC
{
public:
  static InputParameters validParams();

  ADFVMassFractionNuclideInflowBC(const InputParameters & parameters);

protected:
  virtual ADReal computeQpResidual() override;

  const unsigned int _mesh_dims;

  // Density of the bulk fluid.
  const Moose::Functor<ADReal> & _density;

  // The velocity field.
  // X-velocity.
  const Moose::Functor<ADReal> & _vel_u;
  // Y-velocity.
  const Moose::Functor<ADReal> * const _vel_v;
  // Z-velocity.
  const Moose::Functor<ADReal> * const _vel_w;

  const Real _inflow_rate;
}; // class ADFVMassFractionNuclideInflowBC
