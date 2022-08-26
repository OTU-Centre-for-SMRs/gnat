#include "ADIsotopeAdvection.h"

registerMooseObject("GnatApp", ADIsotopeAdvection);

InputParameters
ADIsotopeAdvection::validParams()
{
  auto params = ADIsotopeBase::validParams();
  params.addClassDescription("Computes the advection term for the isotope "
                             "scalar transport equation with a constant "
                             "velocity: "
                             "$(- ( \\vec{\\nabla}\\psi_{j}\\cdot\\vec{v},"
                             "N_{i} )_{V}$. This kernel is stabilized with the "
                             "SUPG method.");

  return params;
}

ADIsotopeAdvection::ADIsotopeAdvection(const InputParameters & parameters)
  : ADIsotopeBase(parameters)
{ }

ADReal
ADIsotopeAdvection::computeQpResidual()
{
  return -1.0 * _grad_test[_i][_qp] * getQpVelocity() * _u[_qp];
}
