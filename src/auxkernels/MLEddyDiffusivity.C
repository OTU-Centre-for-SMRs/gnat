#include "MLEddyDiffusivity.h"

#include "metaphysicl/raw_type.h"

registerMooseObject("GnatApp", MLEddyDiffusivity);

InputParameters
MLEddyDiffusivity::validParams()
{
  InputParameters params = AuxKernel::validParams();
  params.addClassDescription("An auxkernel to compute the eddy diffusivity that appears in the "
                             "RANS form of the trace particle conservation equations.");

  params.addRequiredParam<MooseFunctorName>("u", "The velocity in the x direction.");
  params.addParam<MooseFunctorName>("v", "The velocity in the y direction.");
  params.addParam<MooseFunctorName>("w", "The velocity in the z direction.");

  params.addRequiredParam<MooseFunctorName>("mixing_length", "The turbulent mixing length.");
  params.addParam<Real>("schmidt_number",
                        0.7,
                        "The turbulent Schmidt number that relates the turbulent scalar "
                        "diffusivity to the turbulent momentum diffusivity.");

  // We need some ghosting for the finite volume fields (we use neighbor information to compute
  // gradient)
  params.addParam<unsigned short>("ghost_layers", 2, "The number of layers of elements to ghost.");
  params.addRelationshipManager(
      "ElementSideNeighborLayers",
      Moose::RelationshipManagerType::GEOMETRIC | Moose::RelationshipManagerType::ALGEBRAIC,
      [](const InputParameters & obj_params, InputParameters & rm_params)
      {
        rm_params.set<unsigned short>("layers") = obj_params.get<unsigned short>("ghost_layers");
        rm_params.set<bool>("use_displaced_mesh") = obj_params.get<bool>("use_displaced_mesh");
      });

  return params;
}

MLEddyDiffusivity::MLEddyDiffusivity(const InputParameters & parameters)
  : AuxKernel(parameters),
    _mesh_dims(_subproblem.mesh().dimension()),
    _u(getFunctor<ADReal>("u")),
    _v(isParamValid("v") ? &getFunctor<ADReal>("v") : nullptr),
    _w(isParamValid("w") ? &getFunctor<ADReal>("w") : nullptr),
    _mixing_len(getFunctor<ADReal>("mixing_length")),
    _schmidt_number(getParam<Real>("schmidt_number")),
    _use_qp_arg(dynamic_cast<MooseVariableFE<RealVectorValue> *>(&_var))
{
  if (!_use_qp_arg && !dynamic_cast<MooseVariableFV<Real> *>(&_var))
    paramError(
        "variable",
        "The variable must be a non-vector, non-array finite-volume/finite-element variable.");

  if (isNodal())
    paramError("variable", "This AuxKernel only supports Elemental fields.");

  if (_mesh_dims >= 2 && !_v)
    mooseError(
        "In 2D or 3D the v component of the velocity must be supplied using the 'v' parameter.");
  if (_mesh_dims >= 3 && !_w)
    mooseError("In 3D, the w component of the velocity must be supplied using the 'w' parameter.");
}

Real
MLEddyDiffusivity::computeValue()
{
  constexpr Real offset = 1e-15; // Prevents the explosion of the derivative of sqrt(x).
  Real eddy_diffusivity = 0.0;
  Real symmetric_strain_tensor_norm = 0.0;

  using MetaPhysicL::raw_value;
  if (_use_qp_arg)
  {
    const auto qp_arg = std::make_tuple(_current_elem, _qp, _qrule);

    const auto grad_u = _u.gradient(qp_arg);
    symmetric_strain_tensor_norm += 2.0 * Utility::pow<2>(raw_value(grad_u(0)));
    if (_mesh_dims >= 2)
    {
      const auto grad_v = _v->gradient(qp_arg);
      symmetric_strain_tensor_norm += 2.0 * Utility::pow<2>(raw_value(grad_v(1))) +
                                      Utility::pow<2>(raw_value(grad_v(0)) + raw_value(grad_u(1)));
      if (_mesh_dims >= 3)
      {
        const auto grad_w = _w->gradient(qp_arg);
        symmetric_strain_tensor_norm +=
            2.0 * Utility::pow<2>(raw_value(grad_w(2))) +
            Utility::pow<2>(raw_value(grad_u(2)) + raw_value(grad_w(0))) +
            Utility::pow<2>(raw_value(grad_v(2)) + raw_value(grad_w(1)));
      }
    }

    symmetric_strain_tensor_norm = std::sqrt(symmetric_strain_tensor_norm + offset);

    eddy_diffusivity = symmetric_strain_tensor_norm * raw_value(_mixing_len(qp_arg)) *
                       raw_value(_mixing_len(qp_arg));
  }
  else
  {
    const auto elem_arg = makeElemArg(_current_elem);

    const auto grad_u = _u.gradient(elem_arg);
    symmetric_strain_tensor_norm += 2.0 * Utility::pow<2>(raw_value(grad_u(0)));
    if (_mesh_dims >= 2)
    {
      const auto grad_v = _v->gradient(elem_arg);
      symmetric_strain_tensor_norm += 2.0 * Utility::pow<2>(raw_value(grad_v(1))) +
                                      Utility::pow<2>(raw_value(grad_v(0)) + raw_value(grad_u(1)));
      if (_mesh_dims >= 3)
      {
        const auto grad_w = _w->gradient(elem_arg);
        symmetric_strain_tensor_norm +=
            2.0 * Utility::pow<2>(raw_value(grad_w(2))) +
            Utility::pow<2>(raw_value(grad_u(2)) + raw_value(grad_w(0))) +
            Utility::pow<2>(raw_value(grad_v(2)) + raw_value(grad_w(1)));
      }
    }

    symmetric_strain_tensor_norm = std::sqrt(symmetric_strain_tensor_norm + offset);

    eddy_diffusivity = symmetric_strain_tensor_norm * raw_value(_mixing_len(elem_arg)) *
                       raw_value(_mixing_len(elem_arg));
  }

  return eddy_diffusivity / _schmidt_number;
}
