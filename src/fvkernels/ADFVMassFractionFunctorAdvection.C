#include "ADFVMassFractionFunctorAdvection.h"

#include "NSFVUtils.h"

registerMooseObject("GnatApp", ADFVMassFractionFunctorAdvection);

InputParameters
ADFVMassFractionFunctorAdvection::validParams()
{
  auto params = FVFluxKernel::validParams();
  params.addClassDescription("Computes the advection term for the isotope "
                             "scalar transport equation: "
                             "$(- ( \\vec{\\nabla}\\psi_{j}\\cdot\\vec{v},"
                             "\\rho\\gamma_{i} )_{V}$.");

  params.addRequiredParam<MooseFunctorName>("density",
                                            "The functor name for the density of the bulk fluid.");

  //----------------------------------------------------------------------------
  // Parameters required for advection.
  params.addRequiredParam<MooseFunctorName>(
      "u", "The functor name for the x-component of the velocity.");
  params.addParam<MooseFunctorName>("v", "The functor name for the y-component of the velocity.");
  params.addParam<MooseFunctorName>("w", "The functor name for the z-component of the velocity.");

  params += Moose::FV::interpolationParameters();

  return params;
}

ADFVMassFractionFunctorAdvection::ADFVMassFractionFunctorAdvection(
    const InputParameters & parameters)
  : FVFluxKernel(parameters),
    _mesh_dims(_subproblem.mesh().dimension()),
    _density(getFunctor<ADReal>("density")),
    _vel_u(getFunctor<ADReal>("u")),
    _vel_v(isParamValid("v") ? &getFunctor<ADReal>("v") : nullptr),
    _vel_w(isParamValid("w") ? &getFunctor<ADReal>("w") : nullptr)

{
  if (_mesh_dims >= 2u && !_vel_v)
    mooseError(
        "In 2D or 3D the v component of the velocity must be supplied using the 'v' parameter.");
  if (_mesh_dims >= 3u && !_vel_w)
    mooseError("In 3D, the w component of the velocity must be supplied using the 'w' parameter.");

  const bool need_more_ghosting =
      Moose::FV::setInterpolationMethods(*this, _advected_interp_method, _velocity_interp_method);
  if (need_more_ghosting && _tid == 0)
  {
    adjustRMGhostLayers(std::max((unsigned short)(3), _pars.get<unsigned short>("ghost_layers")));

    // If we need more ghosting, then we are a second-order nonlinear limiting scheme whose stencil
    // is liable to change upon wind-direction change. Consequently we need to tell our problem that
    // it's ok to have new nonzeros which may crop-up after PETSc has shrunk the matrix memory
    getCheckedPointerParam<FEProblemBase *>("_fe_problem_base")
        ->setErrorOnJacobianNonzeroReallocation(false);
  }
}

ADReal
ADFVMassFractionFunctorAdvection::computeQpResidual()
{
  // Get the velocity through central differencing. TODO: More sophisticated interpolation.
  auto face_vel = makeCDFace(*_face_info);
  ADRealVectorValue vel = ADRealVectorValue(_vel_u(face_vel, 0u),
                                            _vel_v ? (*_vel_v)(face_vel, 0u) : 0.0,
                                            _vel_w ? (*_vel_w)(face_vel, 0u) : 0.0);

  // Interpolate the advected quantity to the face.
  const auto face_adv =
      makeFace(*_face_info, Moose::FV::limiterType(_advected_interp_method), vel * _normal > 0.0);
  ADReal u_interface = _var(face_adv, 0u) * _density(face_adv, 0u);

  return _normal * vel * u_interface;
}
