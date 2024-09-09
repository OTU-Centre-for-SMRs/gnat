#include "SASFStreamingJumpIndicator.h"

registerMooseObject("MooseApp", SASFStreamingJumpIndicator);

InputParameters
SASFStreamingJumpIndicator::validParams()
{
  auto params = InternalSideIndicator::validParams();
  params.addClassDescription("An class which uses the jump in particle current as an indicator.");
  params.addRequiredParam<Point>("source_location", "The location of the point source.");

  return params;
}

SASFStreamingJumpIndicator::SASFStreamingJumpIndicator(const InputParameters & parameters)
  : InternalSideIndicator(parameters), _source_location(getParam<Point>("source_location"))
{
}

Real
SASFStreamingJumpIndicator::computeQpIntegral()
{
  auto dir = _q_point[_qp] - _source_location;
  auto mag = std::max(dir.norm(), libMesh::TOLERANCE * libMesh::TOLERANCE);
  const auto omega = dir / mag;

  auto jump = omega * (_grad_u[_qp] - _grad_u_neighbor[_qp]);
  return jump * jump;
}
