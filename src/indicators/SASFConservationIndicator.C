#include "SASFConservationIndicator.h"

registerMooseObject("MooseApp", SASFConservationIndicator);

InputParameters
SASFConservationIndicator::validParams()
{
  auto params = ElementIntegralIndicator::validParams();
  params.addClassDescription(
      "An class which uses the element-integrated lack of local conservation as an indicator.");
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

SASFConservationIndicator::SASFConservationIndicator(const InputParameters & parameters)
  : ElementIntegralIndicator(parameters),
    _source_location(getParam<Point>("source_location")),
    _group_index(getParam<unsigned int>("group_index")),
    _sigma_t_g(getADMaterialProperty<std::vector<Real>>(getParam<std::string>("transport_system") +
                                                        "total_xs_g"))
{
}

Real
SASFConservationIndicator::computeQpIntegral()
{
  auto dir = _q_point[_qp] - _source_location;
  auto mag = std::max(dir.norm(), libMesh::TOLERANCE * libMesh::TOLERANCE);
  const auto omega = dir / mag;

  const auto total = MetaPhysicL::raw_value(_sigma_t_g[_qp][_group_index]);

  auto err = omega * _grad_u[_qp] + total * _u[_qp];

  return err * err;
}
