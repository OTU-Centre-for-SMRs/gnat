#include "ADIsotopeInflowBC.h"

registerMooseObject("GnatApp", ADIsotopeInflowBC);

InputParameters
ADIsotopeInflowBC::validParams()
{
  auto params = ADIsotopeBaseBC::validParams();
  params.addClassDescription("Computes the inflow boundary condition contribution for the isotope "
                             "mass transport equation stabilized with the SUPG method. The weak "
                             "form is given by $(\\phi + \\tau\\vec{v}\\cdot\\vec{\\nabla}\\phi, "
                             "N\\vec{n}\\cdot\\hat{n}\\)$.");
  params.addRequiredRangeCheckedParam<Real>(
      "inflow_rate",
      "inflow_rate >= 0.0",
      "The mass flux or isotope flux entering the domain across the given boundary. The units must "
      "remain consistent with the units of the nonlinear variable.");

  return params;
}

ADIsotopeInflowBC::ADIsotopeInflowBC(const InputParameters & parameters)
  : ADIsotopeBaseBC(parameters), _inflow_rate(getParam<Real>("inflow_rate"))
{
}

ADReal
ADIsotopeInflowBC::computeQpResidual()
{
  return computeQpTests() * getQpVelocity() * _normals[_qp] * _inflow_rate;
}
