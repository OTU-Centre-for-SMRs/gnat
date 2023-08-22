#include "SNRemoval.h"

registerMooseObject("GnatApp", SNRemoval);

InputParameters
SNRemoval::validParams()
{
  auto params = Kernel::validParams();
  params.addClassDescription("Computes the collision term for the "
                             "discrete ordinates neutron transport equation. "
                             "The weak form is given by "
                             "$(\\psi_{j}, \\Sigma_{t,g} \\Psi_{g, n}^{k})$. "
                             "This kernel should not be exposed to the user, "
                             "instead being enabled through a transport action.");
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

SNRemoval::SNRemoval(const InputParameters & parameters)
  : Kernel(parameters),
    _group_index(getParam<unsigned int>("group_index")),
    _sigma_t_g(getADMaterialProperty<std::vector<Real>>(getParam<std::string>("transport_system") +
                                                        "total_xs_g"))
{
}

Real
SNRemoval::computeQpResidual()
{
  return _test[_i][_qp] * MetaPhysicL::raw_value(_sigma_t_g[_qp][_group_index]) * _u[_qp];
}

Real
SNRemoval::computeQpJacobian()
{
  return _test[_i][_qp] * MetaPhysicL::raw_value(_sigma_t_g[_qp][_group_index]) * _phi[_j][_qp];
}
