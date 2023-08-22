#include "INSFVMassFractionScalarFieldAdvection.h"

registerMooseObject("NavierStokesApp", INSFVMassFractionScalarFieldAdvection);

InputParameters
INSFVMassFractionScalarFieldAdvection::validParams()
{
  auto params = INSFVAdvectionKernel::validParams();
  params.addClassDescription("Advects an arbitrary quantity, the associated nonlinear 'variable'.");
  params.addRequiredParam<MooseFunctorName>("density",
                                            "The functor name for the density of the bulk fluid.");
  return params;
}

INSFVMassFractionScalarFieldAdvection::INSFVMassFractionScalarFieldAdvection(
    const InputParameters & params)
  : INSFVAdvectionKernel(params), _density(getFunctor<ADReal>("density"))
{
}

ADReal
INSFVMassFractionScalarFieldAdvection::computeQpResidual()
{
  const auto v = _rc_vel_provider.getVelocity(_velocity_interp_method, *_face_info, 0u, _tid);
  // The mass fraction of the scalar quantity.
  const auto var_face = _var(makeFace(*_face_info,
                                      limiterType(_advected_interp_method),
                                      MetaPhysicL::raw_value(v) * _normal > 0),
                             0u);

  // The density of the bulk fluid.
  const auto density_face = _density(makeFace(*_face_info,
                                              limiterType(_advected_interp_method),
                                              MetaPhysicL::raw_value(v) * _normal > 0),
                                     0u);

  return _normal * v * var_face * density_face;
}
