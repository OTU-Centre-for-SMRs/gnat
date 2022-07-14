#pragma once

#include "AuxKernel.h"

#include "GaussAngularQuadrature.h"
#include "RealSphericalHarmonics.h"

class NeutronFluxMoment : public AuxKernel
{
public:
  static InputParameters validParams();

  NeutronFluxMoment(const InputParameters & parameters);

protected:
  void cartesianToSpherical(const RealVectorValue & ordinate, Real & mu,
                            Real & omega);

  virtual Real computeValue() override;

  std::vector<const ADVariableValue *> _flux_ordinates;

  const unsigned int _degree;
  const int _order;

  const ADMaterialProperty<std::vector<RealVectorValue>> & _quadrature_directions;
  const ADMaterialProperty<std::vector<Real>> & _quadrature_weights;
  const MaterialProperty<MajorAxis> & _axis;
}; // class NeutronFluxMoment
