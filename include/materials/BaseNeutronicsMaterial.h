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

  const unsigned int _harmonics_order;

  GaussAngularQuadrature _quadrature_set;

  std::vector<RealVectorValue> _testVector;
  std::vector<Real> _testVector_2;

  ADMaterialProperty<std::vector<RealVectorValue>> & _quadrature_directions;
  ADMaterialProperty<std::vector<Real>> & _quadrature_weights;
};
