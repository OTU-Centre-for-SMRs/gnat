#include "ADNeutronMatSourceBC.h"

registerMooseObject("GnatApp", ADNeutronMatSourceBC);

InputParameters
ADNeutronMatSourceBC::validParams()
{
  auto params = ADIntegratedBC::validParams();
  params.addClassDescription("Computes the surface source boundary condition "
                             "with the weak form given by "
                             "$\\langle \\psi_{j},\\, \\hat{n}\\cdot\\hat{\\Omega}\\Psi_{inc,\\, g}\\rangle_{\\Gamma_{i}}$, "
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

ADNeutronMatSourceBC::ADNeutronMatSourceBC(const InputParameters & parameters)
  : ADIntegratedBC(parameters)
  , _ordinate_index(getParam<unsigned int>("ordinate_index"))
  , _directions(getADMaterialProperty<std::vector<RealVectorValue>>("directions"))
  , _surface_source(getADMaterialProperty<std::vector<Real>>("surface_source"))
{ }

ADReal
ADNeutronMatSourceBC::computeQpResidual()
{
  if (_ordinate_index >= _directions[_qp].size()
      || _ordinate_index >= _surface_source[_qp].size())
  {
    mooseError("The ordinates index exceeds the number of quadrature points "
               "and/or source values.");
  }

  ADReal res = 0.0;
  ADReal n_dot_omega = _directions[_qp][_ordinate_index] * _normals[_qp];
  if (n_dot_omega > 0.0)
    res += _u[_qp] * n_dot_omega * _test[_i][_qp];
  else
    res += _surface_source[_qp][_ordinate_index] * n_dot_omega * _test[_i][_qp];

  return res;
}
