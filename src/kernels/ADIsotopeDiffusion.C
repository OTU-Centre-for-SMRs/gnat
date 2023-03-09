#include "ADIsotopeDiffusion.h"

registerMooseObject("GnatApp", ADIsotopeDiffusion);

InputParameters
ADIsotopeDiffusion::validParams()
{
  auto params = ADIsotopeBase::validParams();
  params.addClassDescription("Computes the diffusion term for the isotope "
                             "scalar transport equation: "
                             "$( \\vec{\\nabla}\\psi_{j}, "
                             "D_{i}\\vec{\\nabla}N_{i} )_{V}$. This kernel "
                             "expects a diffusion coefficient from the "
                             "material system.");

  return params;
}

ADIsotopeDiffusion::ADIsotopeDiffusion(const InputParameters & parameters)
  : ADIsotopeBase(parameters),
    _grad_mat_diff(getADMaterialProperty<RealVectorValue>(
        "grad_isotope_diff_" + Moose::stringify(getParam<NonlinearVariableName>("variable"))))
{
}

ADReal
ADIsotopeDiffusion::computeQpResidual()
{
  ADReal eddy_diff = 0.0;
  if (_eddy_diffusivity)
    eddy_diff = (*_eddy_diffusivity)[_qp];

  const ADRealVectorValue vel = getQpVelocity();
  // Unstabilized contribution.
  ADReal res = _grad_test[_i][_qp] * (_mat_diff[_qp] + eddy_diff) * _grad_u[_qp];
  // Stabilizing upwind diffusion.
  res -= computeQpTau() * vel * _grad_test[_i][_qp] * (_grad_mat_diff[_qp]) * _grad_u[_qp];

  return res;
}
