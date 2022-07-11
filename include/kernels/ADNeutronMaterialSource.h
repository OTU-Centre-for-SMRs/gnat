#pragma once

#include "ADKernel.h"

#include "GaussAngularQuadrature.h"

class ADNeutronMaterialSource : public ADKernel
{
public:
  static InputParameters validParams();

  ADNeutronMaterialSource(const InputParameters & parameters);

protected:
  virtual ADReal computeQpResidual() override;

  void cartesianToSpherical(const RealVectorValue & ordinate,
                            Real & mu, Real & omega);

  const unsigned int _ordinate_index; // n
  const GaussAngularQuadrature::MajorAxis _axis;

  const ADMaterialProperty<std::vector<RealVectorValue>> & _directions;
  const ADMaterialProperty<std::vector<Real>> & _source_moments;
};
