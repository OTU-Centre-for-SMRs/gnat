#include "SASFAnalyticalFluxBC.h"

registerMooseObject("GnatApp", SASFAnalyticalFluxBC);

InputParameters
SASFAnalyticalFluxBC::validParams()
{
  auto params = IntegratedBC::validParams();
  params.addClassDescription(
      "Imposes a scalar flux boundary condition for the point source uncollided flux "
      "equation stabilized with the self-adjoint flux formalism.");
  params.addRequiredRangeCheckedParam<unsigned int>("group_index",
                                                    "group_index >= 0",
                                                    "The energy group index "
                                                    "of the current angular "
                                                    "flux.");
  params.addRequiredParam<std::vector<Real>>("group_total",
                                             "The macroscopic total "
                                             "cross-sections for all energy "
                                             "groups.");
  params.addRequiredParam<std::vector<Real>>("group_source",
                                             "The group-wise isotropic point source intensity.");
  params.addRequiredParam<Point>("source_location", "The location of the point source.");

  return params;
}

SASFAnalyticalFluxBC::SASFAnalyticalFluxBC(const InputParameters & parameters)
  : IntegratedBC(parameters),
    _group_index(getParam<unsigned int>("group_index")),
    _sigma_t_g(getParam<std::vector<Real>>("group_total")),
    _source_location(getParam<Point>("source_location")),
    _group_source(getParam<std::vector<Real>>("group_source")),
    _dim(_subproblem.mesh().dimension())
{
}

Real
SASFAnalyticalFluxBC::computeQpResidual()
{
  const auto dir = _q_point[_qp] - _source_location;
  const auto r = (_q_point[_qp] - _source_location).norm();
  const auto n_dot_omega = (dir / r) * _normals[_qp];

  Real val = 0.0;
  if (n_dot_omega < 0.0)
  {
    val = _group_source[_group_index] * std::exp(-1.0 * r * _sigma_t_g[_group_index]);
    val /= (_dim == 2u ? (2.0 * libMesh::pi * r) : (4.0 * libMesh::pi * r * r));
  }
  else
    val = _u[_qp];

  return _test[_i][_qp] * n_dot_omega * val;
}
