#include "ADIsotopeDecay.h"

registerMooseObject("GnatApp", ADIsotopeDecay);

InputParameters
ADIsotopeDecay::validParams()
{
  auto params = ADKernel::validParams();
  params.addClassDescription("Computes the radioactive decay sink for the "
                             "isotope scalar transport equation: "
                             "$( \\psi_{j},\\lambda_{i}N_{i} )$.");
  params.addRequiredParam<Real>("decay_const",
                                "The decay constant for this isotope.");

  return params;
}

ADIsotopeDecay::ADIsotopeDecay(const InputParameters & parameters)
  : ADKernel(parameters)
  , _decay_const(getParam<Real>("decay_const"))
{ }

ADReal
ADIsotopeDecay::computeQpResidual()
{
  return _test[_i][_qp] * _decay_const * _u[_qp];
}
