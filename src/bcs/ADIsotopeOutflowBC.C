#include "ADIsotopeOutflowBC.h"

registerMooseObject("GnatApp", ADIsotopeOutflowBC);

InputParameters
ADIsotopeOutflowBC::validParams()
{
  auto params = ADIsotopeBaseBC::validParams();
  params.addClassDescription("Computes the outflow boundary condition contribution for the isotope "
                             "mass transport equation stabilized with the SUPG method. The weak "
                             "form is given by $(\\phi + \\tau\\vec{v}\\cdot\\vec{\\nabla}\\phi, "
                             "N\\vec{n}\\cdot\\hat{n}\\)$.");

  return params;
}

ADIsotopeOutflowBC::ADIsotopeOutflowBC(const InputParameters & parameters)
  : ADIsotopeBaseBC(parameters)
{
}

ADReal
ADIsotopeOutflowBC::computeQpResidual()
{
  ADReal res = 0.0;
  ADRealVectorValue vel = getQpVelocity();
  ADReal n_dot_v = vel * _normals[_qp];
  if (n_dot_v >= 0.0)
    res += n_dot_v * _u[_qp] * _test[_i][_qp];

  return res;
}
