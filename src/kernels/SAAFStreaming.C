#include "SAAFStreaming.h"

registerMooseObject("GnatApp", SAAFStreaming);

InputParameters
SAAFStreaming::validParams()
{
  auto params = SAAFBaseKernel::validParams();
  params.addClassDescription("Computes the streaming term for the "
                             "SAAF discrete ordinates neutron transport equation. "
                             "The weak form is given by "
                             "$(\\nabla\\psi_{j} \\cdot \\hat{\\Omega}, "
                             "\\tau_{g} \\nabla\\Psi_{g, n}^{k} \\cdot "
                             "\\hat{\\Omega} + (\\tau_{g}\\Sigma_{r,\\,g}) - 1)"
                             "\\Psi_{g, n}^{k}$. This kernel should not be "
                             "exposed to the user, instead being enabled "
                             "through a transport action.");
  return params;
}

SAAFStreaming::SAAFStreaming(const InputParameters & parameters) : SAAFBaseKernel(parameters) {}

Real
SAAFStreaming::computeQpResidual()
{
  const auto & omega = _aq.direction(_ordinate_index);

  Real res = MetaPhysicL::raw_value(_saaf_tau[_qp][_group_index]) * omega * _grad_u[_qp] -
             (1.0 - MetaPhysicL::raw_value(_saaf_tau[_qp][_group_index]) *
                        MetaPhysicL::raw_value(_sigma_t_g[_qp][_group_index])) *
                 _u[_qp];
  return _grad_test[_i][_qp] * omega * res;
}

Real
SAAFStreaming::computeQpJacobian()
{
  const auto & omega = _aq.direction(_ordinate_index);

  Real jac = MetaPhysicL::raw_value(_saaf_tau[_qp][_group_index]) * omega * _grad_phi[_j][_qp] -
             (1.0 - MetaPhysicL::raw_value(_saaf_tau[_qp][_group_index]) *
                        MetaPhysicL::raw_value(_sigma_t_g[_qp][_group_index])) *
                 _phi[_j][_qp];
  return _grad_test[_i][_qp] * omega * jac;
}
