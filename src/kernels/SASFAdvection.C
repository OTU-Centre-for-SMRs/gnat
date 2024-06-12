#include "SASFAdvection.h"

registerMooseObject("GnatApp", SASFAdvection);

InputParameters
SASFAdvection::validParams()
{
  auto params = Kernel::validParams();
  params.addClassDescription("The advection contribution for the point source uncollided flux "
                             "equation stabilized with the self-adjoint flux formalism.");
  params.addRequiredRangeCheckedParam<unsigned int>("group_index",
                                                    "group_index >= 0",
                                                    "The energy group index "
                                                    "of the current angular "
                                                    "flux.");
  params.addRequiredParam<Point>("source_location", "The location of the point source.");
  params.addParam<std::string>(
      "transport_system",
      "",
      "Name of the transport system which will consume the provided material properties. If one is "
      "not provided the first transport system will be used.");

  return params;
}

SASFAdvection::SASFAdvection(const InputParameters & parameters)
  : Kernel(parameters),
    _source_location(getParam<Point>("source_location")),
    _div_mult(static_cast<Real>(_subproblem.mesh().dimension() - 1u)),
    _group_index(getParam<unsigned int>("group_index")),
    _sigma_t_g(getADMaterialProperty<std::vector<Real>>(getParam<std::string>("transport_system") +
                                                        "total_xs_g")),
    _saaf_tau(getADMaterialProperty<std::vector<Real>>(getParam<std::string>("transport_system") +
                                                       "saaf_tau"))
{
}

// TODO: Why does the divergence of the direction unit vector cause issues? Negating the divergence
// seems to yield a more correct solution...
Real
SASFAdvection::computeQpResidual()
{
  // SASF.
  auto dir = _q_point[_qp] - _source_location;
  auto mag = std::max(dir.norm(), libMesh::TOLERANCE * libMesh::TOLERANCE);
  const auto omega = dir / mag;
  const auto div_omega = _div_mult / mag;
  const auto tau = MetaPhysicL::raw_value(_saaf_tau[_qp][_group_index]);

  const auto test = omega * _grad_test[_i][_qp];
  const auto grad = tau * omega * _grad_u[_qp];
  const auto lin =
      (tau * (div_omega + MetaPhysicL::raw_value(_sigma_t_g[_qp][_group_index])) - 1.0) * _u[_qp];

  return test * (grad + lin);
}

Real
SASFAdvection::computeQpJacobian()
{
  // SASF.
  auto dir = _q_point[_qp] - _source_location;
  auto mag = std::max(dir.norm(), libMesh::TOLERANCE * libMesh::TOLERANCE);
  const auto omega = dir / mag;
  const auto div_omega = _div_mult / mag;
  const auto tau = MetaPhysicL::raw_value(_saaf_tau[_qp][_group_index]);

  const auto test = omega * _grad_test[_i][_qp];
  const auto grad = tau * omega * _grad_phi[_j][_qp];
  const auto lin =
      (tau * (div_omega + MetaPhysicL::raw_value(_sigma_t_g[_qp][_group_index])) - 1.0) *
      _phi[_j][_qp];

  return test * (grad + lin);
}
