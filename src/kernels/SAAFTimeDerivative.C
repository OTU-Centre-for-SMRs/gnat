#include "SAAFTimeDerivative.h"

// MOOSE includes
#include "Assembly.h"
#include "MooseVariableFE.h"
#include "SystemBase.h"

#include "libmesh/quadrature.h"

registerMooseObject("GnatApp", SAAFTimeDerivative);

InputParameters
SAAFTimeDerivative::validParams()
{
  auto params = SAAFBaseKernel::validParams();
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

SAAFTimeDerivative::SAAFTimeDerivative(const InputParameters & parameters)
  : SAAFBaseKernel(parameters),
    _u_dot(_var.uDot()),
    _du_dot_du(_var.duDotDu()),
    _inv_v_g(getADMaterialProperty<std::vector<Real>>(getParam<std::string>("transport_system") +
                                                      "inv_v_g"))
{
}

Real
SAAFTimeDerivative::computeQpResidual()
{
  return computeQpTests() * MetaPhysicL::raw_value(_inv_v_g[_qp][_group_index]) * _u_dot[_qp];
}

Real
SAAFTimeDerivative::computeQpJacobian()
{
  return computeQpTests() * MetaPhysicL::raw_value(_inv_v_g[_qp][_group_index]) * _du_dot_du[_qp] *
         _phi[_j][_qp];
}
