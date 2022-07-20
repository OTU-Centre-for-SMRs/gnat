#include "ADNeutronStreaming.h"

registerMooseObject("GnatApp", ADNeutronStreaming);

InputParameters
ADNeutronStreaming::validParams()
{
  auto params = ADNeutronBaseKernel::validParams();
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

ADNeutronStreaming::ADNeutronStreaming(const InputParameters & parameters)
  : ADNeutronBaseKernel(parameters)
  , _ordinate_index(getParam<unsigned int>("ordinate_index"))
{ }

ADReal
ADNeutronStreaming::computeQpResidual()
{
  if (_ordinate_index >= _quadrature_set.totalOrder())
    mooseError("The ordinates index exceeds the number of quadrature points.");

  return -1.0 * _grad_test[_i][_qp] * _quadrature_set.direction(_ordinate_index)
         * _u[_qp];
}
