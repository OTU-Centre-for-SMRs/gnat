#include "ADFVNuclideDecaySink.h"

registerMooseObject("GnatApp", ADFVNuclideDecaySink);

InputParameters
ADFVNuclideDecaySink::validParams()
{
  auto params = FVElementalKernel::validParams();
  params.addClassDescription("Computes the radioactive decay sink for the "
                             "isotope scalar transport equation: "
                             "$( \\psi_{j},\\lambda_{i}N_{i} )$.");
  params.addRequiredParam<Real>("decay_const", "The decay constant for this isotope.");

  return params;
}

ADFVNuclideDecaySink::ADFVNuclideDecaySink(const InputParameters & parameters)
  : FVElementalKernel(parameters), _decay_const(getParam<Real>("decay_const"))
{
}

ADReal
ADFVNuclideDecaySink::computeQpResidual()
{
  return _decay_const * _u_functor(makeElemArg(_current_elem), 0u);
}
