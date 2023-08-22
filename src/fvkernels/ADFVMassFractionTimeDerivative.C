#include "ADFVMassFractionTimeDerivative.h"

registerMooseObject("GnatApp", ADFVMassFractionTimeDerivative);

InputParameters
ADFVMassFractionTimeDerivative::validParams()
{
  auto params = FVElementalKernel::validParams();
  params.addClassDescription("Computes the time derivative term for the "
                             "isotope scalar transport equation, where a mass fraction of the bulk "
                             "fluid density is used to represent each nuclide.");
  params.addRequiredParam<MooseFunctorName>("density",
                                            "The functor name for the density of the bulk fluid.");

  return params;
}

ADFVMassFractionTimeDerivative::ADFVMassFractionTimeDerivative(const InputParameters & parameters)
  : FVTimeKernel(parameters), _density(getFunctor<ADReal>("density"))
{
}

ADReal
ADFVMassFractionTimeDerivative::computeQpResidual()
{
  auto elem_args = makeElemArg(_current_elem);

  return _density(elem_args, 0u) * _u_functor.dot(elem_args, 0u) +
         _density.dot(elem_args, 0u) * _u_functor(elem_args, 0u);
}
