#include "ADIsotopeBaseBC.h"

#include "Function.h"

InputParameters
ADIsotopeBaseBC::validParams()
{
  auto params = ADIntegratedBC::validParams();
  params.addClassDescription("A base class which computes the S/U-PG "
                             "stabilization term for the isotope mass "
                             "transport boundary conditions.");

  //----------------------------------------------------------------------------
  // Parameters required for advection and stabilization.
  params.addRequiredParam<MooseFunctorName>(
      "u", "The functor name for the x-component of the velocity.");
  params.addParam<MooseFunctorName>("v", "The functor name for the y-component of the velocity.");
  params.addParam<MooseFunctorName>("w", "The functor name for the z-component of the velocity.");

  params.addRequiredParam<MooseFunctorName>("density",
                                            "The functor name for the density of the bulk fluid.");

  return params;
}

ADIsotopeBaseBC::ADIsotopeBaseBC(const InputParameters & parameters)
  : ADIntegratedBC(parameters),
    _mesh_dims(_subproblem.mesh().dimension()),
    _vel_u(getFunctor<ADReal>("u")),
    _vel_v(isParamValid("v") ? &getFunctor<ADReal>("v") : nullptr),
    _vel_w(isParamValid("w") ? &getFunctor<ADReal>("w") : nullptr),
    _density(getFunctor<ADReal>("density"))
{
  if (_mesh_dims >= 2u && !_vel_v)
    mooseError(
        "In 2D or 3D the v component of the velocity must be supplied using the 'v' parameter.");
  if (_mesh_dims >= 3u && !_vel_w)
    mooseError("In 3D, the w component of the velocity must be supplied using the 'w' parameter.");
}

ADRealVectorValue
ADIsotopeBaseBC::getQpVelocity()
{
  const auto qp_arg = std::make_tuple(_current_elem, _qp, _qrule);
  return ADRealVectorValue(_vel_u(qp_arg, 0u),
                           _vel_v ? (*_vel_v)(qp_arg, 0u) : 0.0,
                           _vel_w ? (*_vel_w)(qp_arg, 0u) : 0.0);
}
