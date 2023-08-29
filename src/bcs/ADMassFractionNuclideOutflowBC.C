#include "ADMassFractionNuclideOutflowBC.h"

registerMooseObject("GnatApp", ADMassFractionNuclideOutflowBC);

InputParameters
ADMassFractionNuclideOutflowBC::validParams()
{
  auto params = ADIsotopeBaseBC::validParams();
  params.addClassDescription("Computes the outflow boundary condition contribution for the isotope "
                             "mass transport equation stabilized with the SUPG method. The weak "
                             "form is given by $(\\phi + \\tau\\vec{v}\\cdot\\vec{\\nabla}\\phi, "
                             "N\\vec{n}\\cdot\\hat{n}\\)$.");

  return params;
}

ADMassFractionNuclideOutflowBC::ADMassFractionNuclideOutflowBC(const InputParameters & parameters)
  : ADIsotopeBaseBC(parameters)
{
}

ADReal
ADMassFractionNuclideOutflowBC::computeQpResidual()
{
  auto qp_arg = Moose::ElemSideQpArg();
  qp_arg.elem = _current_elem;
  qp_arg.side = _current_side;
  qp_arg.qp = _qp;
  qp_arg.point = _q_point[_qp];

  ADReal res = 0.0;
  ADRealVectorValue vel = getQpVelocity();
  ADReal n_dot_v = vel * _normals[_qp];
  if (n_dot_v >= 0.0)
    res += n_dot_v * _u[_qp] * _test[_i][_qp] * _density(qp_arg, 0u);

  return res;
}
