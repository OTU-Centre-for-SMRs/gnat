#include "ADMassFractionNuclideTimeDerivative.h"

registerMooseObject("GnatApp", ADMassFractionNuclideTimeDerivative);

InputParameters
ADMassFractionNuclideTimeDerivative::validParams()
{
  auto params = ADIsotopeBase::validParams();
  params.addClassDescription(
      "Computes the time derivative of the isotope mass transport equation "
      "stabilized with the SUPG method. The weak form is given by $(\\phi + "
      "\\tau\\vec{v}\\cdot\\vec{\\nabla}\\phi,\\, \\frac{\\partial}{\\partial t}N)$.");
  params.set<MultiMooseEnum>("vector_tags") = "time";
  params.set<MultiMooseEnum>("matrix_tags") = "system time";

  return params;
}

ADMassFractionNuclideTimeDerivative::ADMassFractionNuclideTimeDerivative(
    const InputParameters & parameters)
  : ADIsotopeBase(parameters), _u_dot(_var.adUDot())
{
}

ADReal
ADMassFractionNuclideTimeDerivative::computeQpResidual()
{
  auto qp_args = Moose::ElemQpArg();
  qp_args.elem = _current_elem;
  qp_args.qp = _qp;
  qp_args.qrule = _qrule;
  qp_args.point = _q_point[_qp];

  return computeQpTests() *
         (_density(qp_args, 0u) * _u_dot[_qp] + _density.dot(qp_args, 0u) * _u[_qp]);
}
