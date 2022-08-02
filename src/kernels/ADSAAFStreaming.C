#include "ADSAAFStreaming.h"

registerMooseObject("GnatApp", ADSAAFStreaming);

InputParameters
ADSAAFStreaming::validParams()
{
  auto params = ADSAAFBaseKernel::validParams();
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

ADSAAFStreaming::ADSAAFStreaming(const InputParameters & parameters)
  : ADSAAFBaseKernel(parameters)
{ }

ADReal
ADSAAFStreaming::computeQpResidual()
{
  if (_ordinate_index >= _quadrature_set.totalOrder())
    mooseError("The ordinates index exceeds the number of quadrature points.");

  const auto tau = computeQPTau();
  const auto & r_g = _sigma_r_g[_qp][_group_index];
  const auto omega = _quadrature_set.direction(_ordinate_index);

  ADReal res = tau * omega * _grad_u[_qp] - (1.0 - tau * r_g) * _u[_qp];

  return _grad_test[_i][_qp] * omega * res;
}
