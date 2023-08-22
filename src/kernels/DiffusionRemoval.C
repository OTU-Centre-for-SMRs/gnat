#include "DiffusionRemoval.h"

registerMooseObject("GnatApp", DiffusionRemoval);

InputParameters
DiffusionRemoval::validParams()
{
  auto params = Kernel::validParams();
  params.addClassDescription("Computes the collision term for the current group of the neutron "
                             "transport equation simplified with the "
                             "diffusion approximation. The weak form is given by "
                             "$(\\psi_{j}, \\Sigma_{r,g} \\Phi_{g}^{k})$. "
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

DiffusionRemoval::DiffusionRemoval(const InputParameters & parameters)
  : Kernel(parameters),
    _group_index(getParam<unsigned int>("group_index")),
    _sigma_r_g(getADMaterialProperty<std::vector<Real>>(getParam<std::string>("transport_system") +
                                                        "removal_xs_g"))
{
}

Real
DiffusionRemoval::computeQpResidual()
{
  return _test[_i][_qp] * MetaPhysicL::raw_value(_sigma_r_g[_qp][_group_index]) * _u[_qp];
}

Real
DiffusionRemoval::computeQpJacobian()
{
  return _test[_i][_qp] * MetaPhysicL::raw_value(_sigma_r_g[_qp][_group_index]) * _phi[_j][_qp];
}
