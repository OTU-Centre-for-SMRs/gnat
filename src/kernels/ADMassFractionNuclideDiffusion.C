#include "ADMassFractionNuclideDiffusion.h"

registerMooseObject("GnatApp", ADMassFractionNuclideDiffusion);

InputParameters
ADMassFractionNuclideDiffusion::validParams()
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

ADMassFractionNuclideDiffusion::ADMassFractionNuclideDiffusion(const InputParameters & parameters)
  : ADIsotopeBase(parameters),
    _hessian_u(_var.adSecondSln()),
    _mat_diff(getFunctor<ADReal>("isotope_diff_" +
                                 Moose::stringify(getParam<NonlinearVariableName>("variable")))),
    _grad_mat_diff(getFunctor<ADRealVectorValue>(
        "grad_isotope_diff_" + Moose::stringify(getParam<NonlinearVariableName>("variable"))))
{
}

ADReal
ADMassFractionNuclideDiffusion::computeQpResidual()
{
  const auto qp_arg = std::make_tuple(_current_elem, _qp, _qrule);

  // SUPG stabilizing test functions.
  auto supg = _supg_tau(qp_arg, 0u) * getQpVelocity() * _grad_test[_i][_qp];

  // Unstabilized contribution.
  ADReal res = _grad_test[_i][_qp] * _mat_diff(qp_arg, 0u) *
               (_density(qp_arg, 0u) * _grad_u[_qp] + _density.gradient(qp_arg, 0u) * _u[_qp]);

  // Upwind contribution 1.
  res -= supg * _grad_mat_diff(qp_arg, 0u) *
         (_density.gradient(qp_arg, 0u) * _u[_qp] + _density(qp_arg, 0u) * _grad_u[_qp]);

  // Laplacian of the mass fraction.
  ADReal laplacian = _hessian_u[_qp](0u, 0u);
  if (_mesh_dims > 1)
  {
    laplacian += _hessian_u[_qp](1u, 1u);
    if (_mesh_dims > 2)
      laplacian += _hessian_u[_qp](2u, 2u);
  }
  // TODO: Laplacian of the density.
  ADReal density_laplacian = 0.0;

  // Upwind contribution 2. Currently makes the assumption that laplacian(density) = 0.
  res -= supg * _mat_diff(qp_arg, 0u) *
         (density_laplacian + 2.0 * _grad_u[_qp] * _density.gradient(qp_arg, 0u) +
          _density(qp_arg, 0u) * laplacian);

  return res;
}
