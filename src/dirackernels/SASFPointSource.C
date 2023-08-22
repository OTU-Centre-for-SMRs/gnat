#include "SASFPointSource.h"

registerMooseObject("GnatApp", SASFPointSource);

InputParameters
SASFPointSource::validParams()
{
  auto params = DiracKernel::validParams();
  params.addClassDescription("The source contribution for the point source uncollided flux "
                             "equation stabilized with the self-adjoint flux formalism.");
  params.addRequiredParam<Real>("group_source",
                                "The particle source for this spectral energy group.");
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

SASFPointSource::SASFPointSource(const InputParameters & parameters)
  : DiracKernel(parameters),
    _source_location(getParam<Point>("source_location")),
    _group_source(getParam<Real>("group_source")),
    _group_index(getParam<unsigned int>("group_index")),
    _saaf_tau(getADMaterialProperty<std::vector<Real>>(getParam<std::string>("transport_system") +
                                                       "saaf_tau"))
{
}

void
SASFPointSource::addPoints()
{
  addPoint(_source_location);
}

Real
SASFPointSource::computeQpResidual()
{
  // SASF.
  const auto dir = _q_point[_qp] - _source_location;
  const auto mag = std::max(dir.norm(), libMesh::TOLERANCE * libMesh::TOLERANCE);
  const auto omega = dir / mag;

  ///*
  auto test = _test[_i][_qp] +
              (MetaPhysicL::raw_value(_saaf_tau[_qp][_group_index]) * omega * _grad_test[_i][_qp]);
  return -1.0 * test * _group_source;
  //*/
  // return -1.0 * _test[_i][_qp] * _group_source;
}
