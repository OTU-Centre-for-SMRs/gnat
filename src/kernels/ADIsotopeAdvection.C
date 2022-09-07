#include "ADIsotopeAdvection.h"

registerMooseObject("GnatApp", ADIsotopeAdvection);

InputParameters
ADIsotopeAdvection::validParams()
{
  auto params = ADIsotopeBase::validParams();
  params.addClassDescription("Computes the advection term for the isotope "
                             "scalar transport equation: "
                             "$(- ( \\vec{\\nabla}\\psi_{j}\\cdot\\vec{v},"
                             "N_{i} )_{V}$. This kernel is stabilized with the "
                             "SUPG method.");

  return params;
}

ADIsotopeAdvection::ADIsotopeAdvection(const InputParameters & parameters)
  : ADIsotopeBase(parameters)
{
}

ADReal
ADIsotopeAdvection::computeQpResidual()
{
  const ADRealVectorValue vel = getQpVelocity();
  // Unstabilized contribution.
  ADReal res = -1.0 * _grad_test[_i][_qp] * vel * _u[_qp];
  // Stabilizing upwind diffusion.
  res += computeQpTau() * vel * _grad_test[_i][_qp] * vel * _grad_u[_qp];

  return res;
}
