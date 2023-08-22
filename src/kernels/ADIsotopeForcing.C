#include "ADIsotopeForcing.h"

#include "Function.h"

registerMooseObject("GnatApp", ADIsotopeForcing);

InputParameters
ADIsotopeForcing::validParams()
{
  auto params = ADIsotopeBase::validParams();
  params.addClassDescription("A forcing function to test isotope kernels.");

  params.addRequiredParam<FunctionName>("forcing", "The name of the forcing function.");

  return params;
}

ADIsotopeForcing::ADIsotopeForcing(const InputParameters & parameters)
  : ADIsotopeBase(parameters), _forcing(getFunctionByName(getParam<FunctionName>("forcing")))
{
}

ADReal
ADIsotopeForcing::computeQpResidual()
{
  return -1.0 * computeQpTests() * _forcing.value(_t, _q_point[_qp]);
}
