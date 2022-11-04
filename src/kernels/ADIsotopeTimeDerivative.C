#include "ADIsotopeTimeDerivative.h"

registerMooseObject("GnatApp", ADIsotopeTimeDerivative);

InputParameters
ADIsotopeTimeDerivative::validParams()
{
  auto params = ADIsotopeBase::validParams();
  params.addClassDescription(
      "Computes the time derivative of the isotope mass transport equation "
      "stabilized with the SUPG method. The weak form is given by $(\\phi + "
      "\\tau\\vec{v}\\cdot\\vec{\\nabla}\\phi,\\, \\frac{\\partial}{\\partial t}N)$.");
  params.set<MultiMooseEnum>("vector_tags") = "time";
  params.set<MultiMooseEnum>("matrix_tags") = "system time";

  return params;
}

ADIsotopeTimeDerivative::ADIsotopeTimeDerivative(const InputParameters & parameters)
  : ADIsotopeBase(parameters), _u_dot(_var.adUDot())
{
}

ADReal
ADIsotopeTimeDerivative::computeQpResidual()
{
  return computeQpTests() * _u_dot[_qp];
}
