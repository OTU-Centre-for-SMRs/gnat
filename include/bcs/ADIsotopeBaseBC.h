#pragma once

#include "ADIntegratedBC.h"

class ADIsotopeBaseBC : public ADIntegratedBC
{
public:
  static InputParameters validParams();

  ADIsotopeBaseBC(const InputParameters & parameters);

protected:
  // Helper function to fetch the velocity at _qp;
  ADRealVectorValue getQpVelocity();

  const unsigned int _mesh_dims;

  // The velocity field.
  // X-velocity.
  const Moose::Functor<ADReal> & _vel_u;
  // Y-velocity.
  const Moose::Functor<ADReal> * const _vel_v;
  // Z-velocity.
  const Moose::Functor<ADReal> * const _vel_w;

  // Density of the bulk fluid.
  const Moose::Functor<ADReal> & _density;
};
