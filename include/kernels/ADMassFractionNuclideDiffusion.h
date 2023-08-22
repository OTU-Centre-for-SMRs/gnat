#pragma once

#include "ADIsotopeBase.h"

class ADMassFractionNuclideDiffusion : public ADIsotopeBase
{
public:
  static InputParameters validParams();

  ADMassFractionNuclideDiffusion(const InputParameters & parameters);

protected:
  virtual ADReal computeQpResidual() override;

  // The Hessian matrix of _u.
  const ADTemplateVariableSecond<Real> & _hessian_u;

  // Diffusion coefficients for stabilization.
  const Moose::Functor<ADReal> & _mat_diff;
  const Moose::Functor<ADRealVectorValue> & _grad_mat_diff;
}; // class ADMassFractionNuclideDiffusion
