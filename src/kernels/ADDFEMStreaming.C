#include "ADDFEMStreaming.h"

registerMooseObject("GnatApp", ADDFEMStreaming);

InputParameters
ADDFEMStreaming::validParams()
{
  auto params = ADSNBaseKernel::validParams();
  params.addClassDescription("Computes the streaming term for the "
                             "discrete ordinates neutron transport equation. "
                             "The weak form is given by "
                             "$-(\\nabla \\psi_{j}\\cdot\\vec{\\Omega}, "
                             "\\Psi_{g, n}^{k})$. "
                             "This kernel should not be exposed to the user, "
                             "instead being enabled through a transport action.");
  params.addRequiredRangeCheckedParam<unsigned int>("ordinate_index",
                                                    "ordinate_index >= 0",
                                                    "The discrete ordinate index "
                                                    "of the current angular "
                                                    "flux.");

  return params;
}

ADDFEMStreaming::ADDFEMStreaming(const InputParameters & parameters)
  : ADSNBaseKernel(parameters)
  , _ordinate_index(getParam<unsigned int>("ordinate_index"))
{ }

ADReal
ADDFEMStreaming::computeQpResidual()
{
  if (_ordinate_index >= _quadrature_set.totalOrder())
    mooseError("The ordinates index exceeds the number of quadrature points.");

  ADReal res = -1.0 * _grad_test[_i][_qp] * _quadrature_set.direction(_ordinate_index)
               * _u[_qp];

  return res;
}
