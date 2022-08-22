#include "ADIsotopeConstAdvection.h"

registerMooseObject("GnatApp", ADIsotopeConstAdvection);

InputParameters
ADIsotopeConstAdvection::validParams()
{
  auto params = ADKernel::validParams();
  params.addClassDescription("Computes the advection term for the isotope "
                             "scalar transport equation with a constant "
                             "velocity: "
                             "$(- ( \\vec{\\nabla}\\psi_{j}\\cdot\\vec{v},"
                             "N_{i} )_{V}$. This kernel is unstable and can "
                             "only be used with discontinuous finite "
                             "elements.");

  params.addRequiredParam<RealVectorValue>("velocity", "The velocity vector.");

  return params;
}

ADIsotopeConstAdvection::ADIsotopeConstAdvection(const InputParameters & parameters)
  : ADKernel(parameters)
  , _vel(getParam<RealVectorValue>("velocity"))
{ }

ADReal
ADIsotopeConstAdvection::computeQpResidual()
{
  return -1.0 * _grad_test[_i][_qp] * _vel * _u[_qp];
}
