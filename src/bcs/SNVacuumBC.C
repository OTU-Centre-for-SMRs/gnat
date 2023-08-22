#include "SNVacuumBC.h"

registerMooseObject("GnatApp", SNVacuumBC);

InputParameters
SNVacuumBC::validParams()
{
  auto params = SNBaseBC::validParams();
  params.addClassDescription("Computes the vacuum boundary condition with a "
                             "weak form given by "
                             "$\\langle \\psi_{j},\\, 0\\rangle_{\\Gamma_{v}}$, "
                             "$\\hat{n}\\cdot\\hat{\\Omega} \\leq 0$. "
                             "This BC should not be exposed to the user, "
                             "instead being enabled through a transport action.");
  params.addRequiredRangeCheckedParam<unsigned int>("ordinate_index",
                                                    "ordinate_index >= 0",
                                                    "The discrete ordinate index "
                                                    "of the current angular "
                                                    "flux.");

  return params;
}

SNVacuumBC::SNVacuumBC(const InputParameters & parameters)
  : SNBaseBC(parameters), _ordinate_index(getParam<unsigned int>("ordinate_index"))
{
  if (_ordinate_index >= _aq.totalOrder())
    mooseError("The ordinates index exceeds the number of quadrature points.");
}

Real
SNVacuumBC::computeQpResidual()
{
  const auto n_dot_omega = _aq.direction(_ordinate_index) * _normals[_qp];
  return n_dot_omega >= 0.0 ? n_dot_omega * _u[_qp] * _test[_i][_qp] : 0.0;
}

Real
SNVacuumBC::computeQpJacobian()
{
  const auto n_dot_omega = _aq.direction(_ordinate_index) * _normals[_qp];
  return n_dot_omega >= 0.0 ? n_dot_omega * _phi[_j][_qp] * _test[_i][_qp] : 0.0;
}
