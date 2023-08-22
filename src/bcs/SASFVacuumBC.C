#include "SASFVacuumBC.h"

registerMooseObject("GnatApp", SASFVacuumBC);

InputParameters
SASFVacuumBC::validParams()
{
  auto params = IntegratedBC::validParams();
  params.addClassDescription(
      "Imposes a vacuum boundary condition for the point source uncollided flux "
      "equation stabilized with the self-adjoint flux formalism.");
  params.addRequiredParam<Point>("source_location", "The location of the point source.");

  return params;
}

SASFVacuumBC::SASFVacuumBC(const InputParameters & parameters)
  : IntegratedBC(parameters), _source_location(getParam<Point>("source_location"))
{
}

Real
SASFVacuumBC::computeQpResidual()
{
  const auto dir = _q_point[_qp] - _source_location;
  const auto mag = std::max(dir.norm(), libMesh::TOLERANCE * libMesh::TOLERANCE);
  const auto omega = dir / mag;
  const auto n_dot_omega = omega * _normals[_qp];

  if (n_dot_omega >= 0.0)
    return _test[_i][_qp] * n_dot_omega * _u[_qp];
  else
    return 0.0;

  // return _test[_i][_qp] * _u[_qp];
}

Real
SASFVacuumBC::computeQpJacobian()
{
  const auto dir = _q_point[_qp] - _source_location;
  const auto mag = std::max(dir.norm(), libMesh::TOLERANCE * libMesh::TOLERANCE);
  const auto omega = dir / mag;
  const auto n_dot_omega = omega * _normals[_qp];

  if (n_dot_omega >= 0.0)
    return _test[_i][_qp] * n_dot_omega * _phi[_j][_qp];
  else
    return 0.0;

  // return _test[_i][_qp] * _phi[_j][_qp];
}
