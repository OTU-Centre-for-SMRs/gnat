#include "SASFTransverseJumpIndicator.h"

registerMooseObject("MooseApp", SASFTransverseJumpIndicator);

InputParameters
SASFTransverseJumpIndicator::validParams()
{
  auto params = InternalSideIndicator::validParams();
  params.addClassDescription("An class which uses the jump in the transverse derivative of the "
                             "streaming direction as an error indicator.");
  params.addRequiredParam<Point>("source_location", "The location of the point source.");

  return params;
}

SASFTransverseJumpIndicator::SASFTransverseJumpIndicator(const InputParameters & parameters)
  : InternalSideIndicator(parameters), _source_location(getParam<Point>("source_location"))
{
}

Real
SASFTransverseJumpIndicator::computeQpIntegral()
{
  auto dir = _q_point[_qp] - _source_location;
  auto mag = std::max(dir.norm(), libMesh::TOLERANCE * libMesh::TOLERANCE);
  const auto omega = dir / mag;

  const auto transverse = omega.cross(Point(0.0, 0.0, 1.0));

  auto jump = transverse * (_grad_u[_qp] - _grad_u_neighbor[_qp]);
  return jump * jump;
}
