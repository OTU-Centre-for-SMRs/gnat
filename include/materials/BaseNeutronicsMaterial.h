#pragma once

#include "Material.h"

#include "GaussAngularQuadrature.h"

class BaseNeutronicsMaterial : public Material
{
public:
  static InputParameters validParams();

  BaseNeutronicsMaterial(const InputParameters & parameters);

protected:
  virtual void computeQpProperties() override;

  GaussAngularQuadrature _quadrature_set;

  ADMaterialProperty<std::vector<RealVectorValue>> & _quadrature_directions;
  ADMaterialProperty<std::vector<Real>> & _quadrature_weights;
}; // class BaseNeutronicsMaterial
