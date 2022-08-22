#include "ADIsotopeDiffusion.h"

registerMooseObject("GnatApp", ADIsotopeDiffusion);

InputParameters
ADIsotopeDiffusion::validParams()
{
  auto params = ADKernel::validParams();
  params.addClassDescription("Computes the diffusion term for the isotope "
                             "scalar transport equation: "
                             "$( \\vec{\\nabla}\\psi_{j}, "
                             "D_{i}\\vec{\\nabla}N_{i} )_{V}$. This kernel "
                             "expects a diffusion coefficient from the "
                             "material system.");
  params.addRequiredParam<unsigned int>("isotope_id", "The ID of the isotope "
                                        "this kernel applies to.");

  return params;
}

ADIsotopeDiffusion::ADIsotopeDiffusion(const InputParameters & parameters)
  : ADKernel(parameters)
  , _mat_diff(getADMaterialProperty<Real>("isotope_diff" + Moose::stringify(getParam<unsigned int>("isotope_id"))))
{ }

ADReal
ADIsotopeDiffusion::computeQpResidual()
{
  return _grad_test[_i][_qp] * _mat_diff[_qp] * _grad_u[_qp];
}
