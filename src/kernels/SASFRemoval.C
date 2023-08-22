#include "SASFRemoval.h"

registerMooseObject("GnatApp", SASFRemoval);

InputParameters
SASFRemoval::validParams()
{
  auto params = Kernel::validParams();
  params.addClassDescription("Computes the total interaction component of the uncollided point "
                             "source equation with the SASF formalism.");
  params.addRequiredRangeCheckedParam<unsigned int>("group_index",
                                                    "group_index >= 0",
                                                    "The energy group index "
                                                    "of the current angular "
                                                    "flux.");
  params.addParam<std::string>(
      "transport_system",
      "",
      "Name of the transport system which will consume the provided material properties. If one is "
      "not provided the first transport system will be used.");

  return params;
}

SASFRemoval::SASFRemoval(const InputParameters & parameters)
  : Kernel(parameters),
    _group_index(getParam<unsigned int>("group_index")),
    _sigma_t_g(getADMaterialProperty<std::vector<Real>>(getParam<std::string>("transport_system") +
                                                        "total_xs_g"))
{
}

Real
SASFRemoval::computeQpResidual()
{
  // SASF.
  return _test[_i][_qp] * MetaPhysicL::raw_value(_sigma_t_g[_qp][_group_index]) * _u[_qp];
}
Real
SASFRemoval::computeQpJacobian()
{
  // SASF.
  return _test[_i][_qp] * MetaPhysicL::raw_value(_sigma_t_g[_qp][_group_index]) * _phi[_j][_qp];
}
