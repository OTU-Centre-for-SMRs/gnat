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
  return computeQpTests() * _decay_const * _u[_qp] *
         _density(std::make_tuple(_current_elem, _qp, _qrule), 0u);
}
