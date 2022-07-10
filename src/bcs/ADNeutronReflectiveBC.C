#include "ADNeutronReflectiveBC.h"

registerMooseObject("GnatApp", ADNeutronReflectiveBC);

InputParameters
ADNeutronReflectiveBC::validParams()
{
  auto params = ADIntegratedBC::validParams();
  params.addClassDescription("Computes the reflected boundary condition with a "
                             "weak form given by "
                             "$\\langle \\psi_{j},\\, \\hat{n}\\cdot\\hat{\\Omega}\\Psi_{r,\\, g}\\rangle_{\\Gamma_{i}}$, "
                             "$\\hat{n}\\cdot\\hat{\\Omega} \\leq 0$. "
                             "This kernel should not be exposed to the user, "
                             "instead being enabled through a transport action.");
  params.addRequiredRangeCheckedParam<unsigned int>("ordinate_index",
                                                    "ordinate_index >= 0",
                                                    "The discrete ordinate index "
                                                    "of the current angular "
                                                    "flux.");
  params.addRequiredCoupledVar("psi_ref", "The variable for "
                               "the reflected angular flux.");

  return params;
}

ADNeutronReflectiveBC::ADNeutronReflectiveBC(const InputParameters & parameters)
  : ADIntegratedBC(parameters)
  , _ordinate_index(getParam<unsigned int>("ordinate_index"))
  , _u_ref(coupledValue("psi_ref"))
  , _albedo(getADMaterialProperty<Real>("boundary_albedo"))
  , _directions(getADMaterialProperty<std::vector<RealVectorValue>>("directions"))
{ }

ADReal
ADNeutronReflectiveBC::computeQpResidual()
{
  if (_ordinate_index >= _directions[_qp].size())
    mooseError("The ordinates index exceeds the number of quadrature points.");

  ADReal res = 0.0;
  ADReal n_dot_omega = _directions[_qp][_ordinate_index] * _normals[_qp];
  if (n_dot_omega > 0.0)
    res += _u[_qp] * n_dot_omega * _test[_i][_qp];
  else
    res += _albedo[_qp] * _u_ref[_qp] * n_dot_omega * _test[_i][_qp];

  return res;
}
