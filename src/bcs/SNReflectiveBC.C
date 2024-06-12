#include "SNReflectiveBC.h"

#include <algorithm>

registerMooseObject("GnatApp", SNReflectiveBC);

bool
SNReflectiveBC::vecEquals(const RealVectorValue & first,
                          const RealVectorValue & second,
                          const Real & tol)
{
  bool comp_1 = std::abs(first(0) - second(0)) <= tol;
  bool comp_2 = std::abs(first(1) - second(1)) <= tol;
  bool comp_3 = std::abs(first(2) - second(2)) <= tol;

  return comp_1 && comp_2 && comp_3;
}

InputParameters
SNReflectiveBC::validParams()
{
  auto params = SNBaseBC::validParams();
  params.addClassDescription(
      "Computes the reflected boundary condition with a "
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
  params.addRequiredCoupledVar("psi_ref", "The variable for the reflected angular fluxes.");
  params.addRequiredParam<std::vector<RealVectorValue>>(
      "unique_normals",
      "The unique normals this boundary has. All of these normals must be identical to one "
      "cartesian axis.");

  return params;
}

SNReflectiveBC::SNReflectiveBC(const InputParameters & parameters)
  : SNBaseBC(parameters),
    _ordinate_index(getParam<unsigned int>("ordinate_index")),
    _unique_normals(getParam<std::vector<RealVectorValue>>("unique_normals"))
{
  if (coupledComponents("psi_ref") != _unique_normals.size())
    mooseError("The number of provided reflected directions and unique normals must match.");

  _reflected_ordinates.reserve(_unique_normals.size());
  for (unsigned int i = 0u; i < _unique_normals.size(); ++i)
  {
    _reflected_ordinates.emplace_back(&coupledValue("psi_ref", i));
    _jvar_map.emplace(coupled("psi_ref", i));
  }
}

unsigned int
SNReflectiveBC::getReflectedIndex(const RealVectorValue & normal)
{
  int refl = -1;
  for (unsigned int i = 0u; i < _unique_normals.size(); ++i)
  {
    if (vecEquals(normal, _unique_normals[i]))
    {
      refl = static_cast<int>(i);
      break;
    }
  }
  if (refl < 0)
  {
    std::string options = "";
    for (const auto & norm : _unique_normals)
      options += Moose::stringify(static_cast<Point>(norm)) + " ";

    mooseError("The normal " + Moose::stringify(static_cast<Point>(normal)) +
               " does not have an associated reflected direction. Options are: " + options);
  }

  return static_cast<unsigned int>(refl);
}

Real
SNReflectiveBC::computeQpResidual()
{
  auto refl = getReflectedIndex(_normals[_qp]);

  Real res = 0.0;
  Real n_dot_omega = _aq.direction(_ordinate_index) * _normals[_qp];
  if (n_dot_omega > 0.0)
    res += _u[_qp] * n_dot_omega * _test[_i][_qp];
  else
    res += (*_reflected_ordinates[refl])[_qp] * n_dot_omega * _test[_i][_qp];

  return res;
}

Real
SNReflectiveBC::computeQpJacobian()
{
  const auto n_dot_omega = _aq.direction(_ordinate_index) * _normals[_qp];
  return n_dot_omega >= 0.0 ? n_dot_omega * _phi[_j][_qp] * _test[_i][_qp] : 0.0;
}

Real
SNReflectiveBC::computeQpOffDiagJacobian(unsigned int jvar)
{
  if (_jvar_map.count(jvar) > 0u)
  {
    const auto n_dot_omega = _aq.direction(_ordinate_index) * _normals[_qp];
    return n_dot_omega < 0.0 ? n_dot_omega * _phi[_j][_qp] * _test[_i][_qp] : 0.0;
  }
  else
    return 0.0;
}
