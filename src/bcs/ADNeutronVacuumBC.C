#include "ADNeutronVacuumBC.h"

registerMooseObject("GnatApp", ADNeutronVacuumBC);

InputParameters
ADNeutronVacuumBC::validParams()
{
  auto params = ADNeutronBaseBC::validParams();
  params.addClassDescription("Computes the vacuum boundary condition with a "
                             "weak form given by "
                             "$\\langle \\psi_{j},\\, 0\\rangle_{\\Gamma_{v}}$, "
                             "$\\hat{n}\\cdot\\hat{\\Omega} \\leq 0$. "
                             "This kernel should not be exposed to the user, "
                             "instead being enabled through a transport action.");
  params.addRequiredRangeCheckedParam<unsigned int>("ordinate_index",
                                                    "ordinate_index >= 0",
                                                    "The discrete ordinate index "
                                                    "of the current angular "
                                                    "flux.");

  return params;
}

ADNeutronVacuumBC::ADNeutronVacuumBC(const InputParameters & parameters)
  : ADNeutronBaseBC(parameters)
  , _ordinate_index(getParam<unsigned int>("ordinate_index"))
{ }

ADReal
ADNeutronVacuumBC::computeQpResidual()
{
  if (_ordinate_index >= _quadrature_set.totalOrder())
    mooseError("The ordinates index exceeds the number of quadrature points.");

  ADReal res = 0.0;
  ADReal n_dot_omega = _quadrature_set.direction(_ordinate_index) * _normals[_qp];
  if (n_dot_omega > 0.0)
    res += _u[_qp] * n_dot_omega * _test[_i][_qp];

  return res;
}
