#pragma once

#include "AuxKernel.h"

// An auxkernel to compute the eddy diffusivity of trace species. Adapted from
// "INSFVMixingLengthScalarDiffusion" to support finite element trace species simulations.
class MLEddyDiffusivity : public AuxKernel
{
public:
  static InputParameters validParams();

  MLEddyDiffusivity(const InputParameters & parameters);

protected:
  virtual Real computeValue() override;

  const unsigned int _mesh_dims;

  // X-velocity.
  const Moose::Functor<ADReal> & _u;
  // Y-velocity.
  const Moose::Functor<ADReal> * const _v;
  // Z-velocity.
  const Moose::Functor<ADReal> * const _w;

  // Turbulent eddy mixing length.
  const Moose::Functor<ADReal> & _mixing_len;

  // Turbulent Schmidt number.
  const Real & _schmidt_number;

  // Whether to use a quadrature-based functor argument, appropriate for finite element
  // evaluations. If false, use a cell-center functor argument appropriate for finite volume
  // calculations
  const bool _use_qp_arg;
}; // class MLEddyDiffusivity
