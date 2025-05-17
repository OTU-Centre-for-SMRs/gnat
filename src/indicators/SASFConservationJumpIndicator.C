#include "SASFConservationJumpIndicator.h"

#include "FEProblemBase.h"

registerMooseObject("MooseApp", SASFConservationJumpIndicator);

InputParameters
SASFConservationJumpIndicator::validParams()
{
  auto params = InternalSideIndicator::validParams();
  params.addClassDescription("An class which uses the jump in particle current as an indicator.");
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

  params.set<Moose::MaterialDataType>("_material_data_type") = Moose::BLOCK_MATERIAL_DATA;

  return params;
}

SASFConservationJumpIndicator::SASFConservationJumpIndicator(const InputParameters & parameters)
  : InternalSideIndicator(parameters),
    _neighbor_material_data(_mi_feproblem.getMaterialData(Moose::NEIGHBOR_MATERIAL_DATA,
                                                          _mi_params.get<THREAD_ID>("_tid"))),
    _source_location(getParam<Point>("source_location")),
    _group_index(getParam<unsigned int>("group_index")),
    _sigma_t_c_g(getGenericMaterialPropertyByName<std::vector<Real>, true>(
        getParam<std::string>("transport_system") + "total_xs_g", _material_data, 0)),
    _sigma_t_n_g(getGenericMaterialPropertyByName<std::vector<Real>, true>(
        getParam<std::string>("transport_system") + "total_xs_g", _neighbor_material_data, 0))
{
}

Real
SASFConservationJumpIndicator::computeQpIntegral()
{
  auto dir = _q_point[_qp] - _source_location;
  auto mag = std::max(dir.norm(), libMesh::TOLERANCE * libMesh::TOLERANCE);
  const auto omega = dir / mag;

  const auto total_c = MetaPhysicL::raw_value(_sigma_t_c_g[_qp][_group_index]);
  auto err_c = omega * _grad_u[_qp] + total_c * _u[_qp];

  const auto total_n = MetaPhysicL::raw_value(_sigma_t_n_g[_qp][_group_index]);
  auto err_n = omega * _grad_u_neighbor[_qp] + total_n * _u_neighbor[_qp];

  auto jump = err_c - err_n;
  return jump * jump * mag * mag * mag * mag;
}
