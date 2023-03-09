#include "ADSAAFTimeDerivative.h"

// MOOSE includes
#include "Assembly.h"
#include "MooseVariableFE.h"
#include "SystemBase.h"

#include "libmesh/quadrature.h"

registerMooseObject("GnatApp", ADSAAFTimeDerivative);

InputParameters
ADSAAFTimeDerivative::validParams()
{
  auto params = ADSAAFBaseKernel::validParams();
  params.addClassDescription("Computes the time derivative term for the "
                             "SAAF discrete ordinates neutron transport equation. "
                             "The weak form is given by $(\\psi_{j} + "
                             "\\tau_{g}\\vec{\\nabla}\\psi_{j}\\cdot"
                             "\\hat{\\Omega}, \\frac{1}{v_{g}}"
                             "\\frac{\\partial}{\\partial t}\\Psi_{g, n}^{k})$. "
                             "This kernel should not be exposed to the user, "
                             "instead being enabled through a transport action.");
  params.set<MultiMooseEnum>("vector_tags") = "time";
  params.set<MultiMooseEnum>("matrix_tags") = "system time";

  return params;
}

ADSAAFTimeDerivative::ADSAAFTimeDerivative(const InputParameters & parameters)
  : ADSAAFBaseKernel(parameters),
    _u_dot(_var.adUDot()),
    _inv_v_g(getADMaterialProperty<std::vector<Real>>(getParam<std::string>("transport_system") + "inv_v_g"))
{
}

ADReal
ADSAAFTimeDerivative::computeQpResidual()
{
  if (_group_index >= _inv_v_g[_qp].size())
    mooseError("The group index exceeds the number of provided neutron speeds.");

  return computeQPTests() * _inv_v_g[_qp][_group_index] * _u_dot[_qp];
}
