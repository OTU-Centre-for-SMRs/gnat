#pragma once

#include "ADKernel.h"

class ADIsotopeBase : public ADKernel
{
public:
  static InputParameters validParams();

  ADIsotopeBase(const InputParameters & parameters);

protected:
  // Helper function to fetch the velocity at a given qp;
  ADRealVectorValue getQpVelocity();

  // Functions for computing SUPG stabilization parameters.
  ADReal computeQpTests();

  const unsigned int _mesh_dims;

  // The velocity field.
  // X-velocity.
  const Moose::Functor<ADReal> & _vel_u;
  // Y-velocity.
  const Moose::Functor<ADReal> * const _vel_v;
  // Z-velocity.
  const Moose::Functor<ADReal> * const _vel_w;

  // The stabilization parameter for the SUPG finite element weak form.
  const Moose::Functor<ADReal> & _supg_tau;

  // Density of the bulk fluid.
  const Moose::Functor<ADReal> & _density;
}; // class ADIsotopeBase
