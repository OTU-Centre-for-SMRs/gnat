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
  const auto qp_arg = std::make_tuple(_current_elem, _qp, _qrule);
  ADReal res = 0.0;
  ADRealVectorValue vel = getQpVelocity();
  ADReal n_dot_v = vel * _normals[_qp];
  if (n_dot_v >= 0.0)
    res += n_dot_v * _u[_qp] * _test[_i][_qp] * _density(qp_arg, 0u);

  return res;
}
