#include "VariableScaleAux.h"

registerMooseObject("GnatApp", VariableScaleAux);

InputParameters
VariableScaleAux::validParams()
{
  auto params = AuxKernel::validParams();
  params.addClassDescription("An auxkernel which scales a field by a scalar value.");
  params.addRequiredParam<bool>("is_fe", "Whether the mass fraction is finite element or not.");
  params.addCoupledVar(
      "unscaled_var",
      "The unscaled variable this kernel should multiply by 'scale_factor' (in variable form).");
  params.addParam<MooseFunctorName>(
      "unscaled_fun",
      "The unscaled variable this kernel should multiply by 'scale_factor' (in functor form).");
  params.addParam<Real>("scale_factor", "The value this kernel should multiply 'unscaled' by.");

  return params;
}

VariableScaleAux::VariableScaleAux(const InputParameters & parameters)
  : AuxKernel(parameters),
    _is_fe(getParam<bool>("is_fe")),
    _to_scale_var(nullptr),
    _to_scale_fun(nullptr),
    _scale_factor(getParam<Real>("scale_factor"))
{
  if (_is_fe)
    _to_scale_var = &adCoupledValue("unscaled_var");
  else
    _to_scale_fun = &getFunctor<ADReal>("unscaled_fun");
}

Real
VariableScaleAux::computeValue()
{
  if (_is_fe)
    return _scale_factor * MetaPhysicL::raw_value((*_to_scale_var)[_qp]);
  else
    return _scale_factor * MetaPhysicL::raw_value((*_to_scale_fun)(makeElemArg(_current_elem), 0u));
}
