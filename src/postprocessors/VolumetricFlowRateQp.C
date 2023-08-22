#include "VolumetricFlowRateQp.h"

registerMooseObject("GnatApp", VolumetricFlowRateQp);

InputParameters
VolumetricFlowRateQp::validParams()
{
  auto params = SideIntegralPostprocessor::validParams();
  params.addClassDescription(
      "Computes the volumetric flow rate of an advected quantity through a sideset.");

  params.addRequiredParam<MooseFunctorName>("advected_quantity",
                                            "The functor name for the advected quantity.");

  params.addRequiredParam<MooseFunctorName>(
      "u", "The functor name for the x-component of the velocity.");
  params.addParam<MooseFunctorName>("v", "The functor name for the y-component of the velocity.");
  params.addParam<MooseFunctorName>("w", "The functor name for the z-component of the velocity.");

  return params;
}

VolumetricFlowRateQp::VolumetricFlowRateQp(const InputParameters & parameters)
  : SideIntegralPostprocessor(parameters),
    _mesh_dims(_subproblem.mesh().dimension()),
    _vel_u(getFunctor<ADReal>("u")),
    _vel_v(isParamValid("v") ? &getFunctor<ADReal>("v") : nullptr),
    _vel_w(isParamValid("w") ? &getFunctor<ADReal>("w") : nullptr),
    _adv_quant(getFunctor<ADReal>("advected_quantity"))
{
  if (_mesh_dims >= 2u && !_vel_v)
    mooseError(
        "In 2D or 3D the v component of the velocity must be supplied using the 'v' parameter.");
  if (_mesh_dims >= 3u && !_vel_w)
    mooseError("In 3D, the w component of the velocity must be supplied using the 'w' parameter.");

  _qp_integration = true;
}

Real
VolumetricFlowRateQp::computeQpIntegral()
{
  const auto qp_arg = std::make_tuple(_current_elem, _qp, _qrule);
  const auto vel = MetaPhysicL::raw_value(ADRealVectorValue(_vel_u(qp_arg, 0u),
                                                            _vel_v ? (*_vel_v)(qp_arg, 0u) : 0.0,
                                                            _vel_w ? (*_vel_w)(qp_arg, 0u) : 0.0));

  return MetaPhysicL::raw_value(_adv_quant(qp_arg, 0u)) * vel * _normals[_qp];
}
