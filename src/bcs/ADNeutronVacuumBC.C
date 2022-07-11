#include "ADNeutronVacuumBC.h"

registerMooseObject("GnatApp", ADNeutronVacuumBC);

InputParameters
ADNeutronVacuumBC::validParams()
{
  auto params = ADIntegratedBC::validParams();
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
  : ADIntegratedBC(parameters)
  , _ordinate_index(getParam<unsigned int>("ordinate_index"))
  , _directions(getADMaterialProperty<std::vector<RealVectorValue>>("directions"))
{ }

ADReal
ADNeutronVacuumBC::computeQpResidual()
{
  if (_ordinate_index >= _directions[_qp].size())
    mooseError("The ordinates index exceeds the number of quadrature points.");

  ADReal res = 0.0;
  ADReal n_dot_omega = _directions[_qp][_ordinate_index] * _normals[_qp];
  if (n_dot_omega > 0.0)
    res += _u[_qp] * n_dot_omega * _test[_i][_qp];

  return res;
}
