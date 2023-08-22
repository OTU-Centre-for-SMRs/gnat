#include "SpecificActivity.h"

registerMooseObject("GnatApp", SpecificActivity);

InputParameters
SpecificActivity::validParams()
{
  auto params = AuxKernel::validParams();
  params.addClassDescription("An auxkernel for the purpose of computing the specific activity of "
                             "a radionuclide concentration for postprocessing.");
  params.addRequiredParam<bool>("is_fe", "Whether the mass fraction is finite element or not.");
  params.addCoupledVar("nuclide_var", "The radionuclide concentration.");
  params.addParam<MooseFunctorName>("nuclide_fun", "The radionuclide concentration.");
  params.addRequiredParam<Real>("decay_const",
                                "The decay constant of the radionuclide scalar field.");

  return params;
}

SpecificActivity::SpecificActivity(const InputParameters & parameters)
  : AuxKernel(parameters),
    _is_fe(getParam<bool>("is_fe")),
    _concentration_var(nullptr),
    _concentration_fun(nullptr),
    _decay_const(getParam<Real>("decay_const"))
{
  if (_is_fe)
    _concentration_var = &adCoupledValue("nuclide_var");
  else
    _concentration_fun = &getFunctor<ADReal>("nuclide_fun");
}

Real
SpecificActivity::computeValue()
{
  if (_is_fe)
    return _decay_const * MetaPhysicL::raw_value((*_concentration_var)[_qp]);
  else
    return _decay_const *
           MetaPhysicL::raw_value((*_concentration_fun)(makeElemArg(_current_elem), 0u));
}
