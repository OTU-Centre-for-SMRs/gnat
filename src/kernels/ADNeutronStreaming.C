#include "ADNeutronStreaming.h"

registerMooseObject("GnatApp", ADNeutronStreaming);

InputParameters
ADNeutronStreaming::validParams()
{
  auto params = ADKernel::validParams();
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
  : ADKernel(parameters)
  , _ordinate_index(getParam<unsigned int>("ordinate_index"))
  , _directions(getMaterialProperty<std::vector<RealVectorValue>>("directions"))
{ }

ADReal
ADNeutronStreaming::computeQpResidual()
{
  if (_ordinate_index >= _directions[_qp].size())
  {
    mooseWarning(Moose::stringify(_directions[_qp].size()));
    mooseError("The ordinates index exceeds the number of quadrature points.");
  }

  return -1.0 * _grad_test[_i][_qp] * _directions[_qp][_ordinate_index] * _u[_qp];
}
