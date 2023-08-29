#include "ADMassFractionNuclideDecaySink.h"

registerMooseObject("GnatApp", ADMassFractionNuclideDecaySink);

InputParameters
ADMassFractionNuclideDecaySink::validParams()
{
  auto params = ADIsotopeBase::validParams();
  params.addClassDescription("Computes the radioactive decay sink for the "
                             "isotope scalar transport equation: "
                             "$( \\psi_{j},\\lambda_{i}N_{i} )$.");
  params.addRequiredParam<Real>("decay_const", "The decay constant for this isotope.");

  return params;
}

ADMassFractionNuclideDecaySink::ADMassFractionNuclideDecaySink(const InputParameters & parameters)
  : ADIsotopeBase(parameters), _decay_const(getParam<Real>("decay_const"))
{
}

ADReal
ADMassFractionNuclideDecaySink::computeQpResidual()
{
  auto qp_args = Moose::ElemQpArg();
  qp_args.elem = _current_elem;
  qp_args.qp = _qp;
  qp_args.qrule = _qrule;
  qp_args.point = _q_point[_qp];

  return computeQpTests() * _decay_const * _u[_qp] * _density(qp_args, 0u);
}
