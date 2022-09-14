#pragma once

#include "Material.h"

class AutoIsotopeMaterial : public Material
{
public:
  static InputParameters validParams();

  AutoIsotopeMaterial(const InputParameters & parameters);

protected:
  virtual void computeQpProperties() override;

  // TODO: Have the diffusion coefficient rely on fluid properties.
  enum class CoefficientType
  {
    Constant = 0u
  } _diff_type;

  const Real _diffusion_coefficient;

  ADMaterialProperty<Real> & _mat_diff;
  ADMaterialProperty<RealVectorValue> & _grad_mat_diff;
};
