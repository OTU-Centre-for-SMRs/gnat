#include "ADSNMatSourceBC.h"

registerMooseObject("GnatApp", ADSNMatSourceBC);

InputParameters
ADSNMatSourceBC::validParams()
{
  auto params = ADSNBaseBC::validParams();
  params.addClassDescription("Computes the surface source boundary condition "
                             "with the weak form given by "
                             "$\\langle \\psi_{j},\\, \\hat{n}\\cdot"
                             "\\hat{\\Omega}\\Psi_{inc,\\, g}\\rangle_{\\Gamma_{i}}$, "
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

ADSNMatSourceBC::ADSNMatSourceBC(const InputParameters & parameters)
  : ADSNBaseBC(parameters)
  , _ordinate_index(getParam<unsigned int>("ordinate_index"))
  , _surface_source(getADMaterialProperty<std::vector<Real>>("surface_source"))
{ }

ADReal
ADSNMatSourceBC::computeQpResidual()
{
  if (_ordinate_index >= _quadrature_set.totalOrder())
    mooseError("The ordinates index exceeds the number of quadrature points.");

  ADReal res = 0.0;
  ADReal n_dot_omega = _quadrature_set.direction(_ordinate_index) * _normals[_qp];
  if (n_dot_omega > 0.0)
    res += _test[_i][_qp] * _u[_qp] * n_dot_omega;
  else
  {
    res += _test[_i][_qp] * _surface_source[_qp][_ordinate_index] * n_dot_omega
           * _symmetry_factor;
  }

  return res;
}
