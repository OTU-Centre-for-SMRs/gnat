#include "ADFVMassFractionNuclideDecaySink.h"

registerMooseObject("GnatApp", ADFVMassFractionNuclideDecaySink);

InputParameters
ADFVMassFractionNuclideDecaySink::validParams()
{
  auto params = FVElementalKernel::validParams();
  params.addClassDescription("Computes the radioactive decay sink for the "
                             "isotope scalar transport equation: "
                             "$( \\psi_{j},\\lambda_{i}N_{i} )$.");
  params.addRequiredParam<MooseFunctorName>("density",
                                            "The functor name for the density of the bulk fluid.");
  params.addRequiredParam<Real>("decay_const", "The decay constant for this isotope.");

  return params;
}

ADFVMassFractionNuclideDecaySink::ADFVMassFractionNuclideDecaySink(
    const InputParameters & parameters)
  : FVElementalKernel(parameters),
    _density(getFunctor<ADReal>("density")),
    _decay_const(getParam<Real>("decay_const"))
{
}

ADReal
ADFVMassFractionNuclideDecaySink::computeQpResidual()
{
  return _decay_const * _u_functor(makeElemArg(_current_elem), 0u) *
         _density(makeElemArg(_current_elem), 0u);
}
