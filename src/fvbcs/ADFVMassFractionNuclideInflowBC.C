#include "ADFVMassFractionNuclideInflowBC.h"

registerMooseObject("GnatApp", ADFVMassFractionNuclideInflowBC);

InputParameters
ADFVMassFractionNuclideInflowBC::validParams()
{
  auto params = FVFluxBC::validParams();
  params.addClassDescription("Computes the inflow boundary condition contribution for the isotope "
                             "mass transport equation stabilized with the SUPG method. The weak "
                             "form is given by $(N\\vec{n}\\cdot\\hat{n}\\)$.");
  params.addRangeCheckedParam<Real>(
      "inflow_rate",
      0.0,
      "inflow_rate >= 0.0",
      "The mass flux or isotope flux entering the domain across the given boundary. The units must "
      "remain consistent with the units of the nonlinear variable.");

  params.addRequiredParam<MooseFunctorName>("density",
                                            "The functor name for the density of the bulk fluid.");

  //----------------------------------------------------------------------------
  // Parameters required for advection.
  params.addRequiredParam<MooseFunctorName>(
      "u", "The functor name for the x-component of the velocity.");
  params.addParam<MooseFunctorName>("v", "The functor name for the y-component of the velocity.");
  params.addParam<MooseFunctorName>("w", "The functor name for the z-component of the velocity.");

  return params;
}

ADFVMassFractionNuclideInflowBC::ADFVMassFractionNuclideInflowBC(const InputParameters & parameters)
  : FVFluxBC(parameters),
    _mesh_dims(_subproblem.mesh().dimension()),
    _density(getFunctor<ADReal>("density")),
    _vel_u(getFunctor<ADReal>("u")),
    _vel_v(isParamValid("v") ? &getFunctor<ADReal>("v") : nullptr),
    _vel_w(isParamValid("w") ? &getFunctor<ADReal>("w") : nullptr),
    _inflow_rate(getParam<Real>("inflow_rate"))
{
  if (_mesh_dims >= 2u && !_vel_v)
    mooseError(
        "In 2D or 3D the v component of the velocity must be supplied using the 'v' parameter.");
  if (_mesh_dims >= 3u && !_vel_w)
    mooseError("In 3D, the w component of the velocity must be supplied using the 'w' parameter.");
}

ADReal
ADFVMassFractionNuclideInflowBC::computeQpResidual()
{
  const auto face = makeCDFace(*_face_info, false);
  const auto vel = ADRealVectorValue(
      _vel_u(face, 0u), _vel_v ? (*_vel_v)(face, 0u) : 0.0, _vel_w ? (*_vel_w)(face, 0u) : 0.0);
  return _normal * vel * _inflow_rate * _density(face, 0u);
}
