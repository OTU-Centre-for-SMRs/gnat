#pragma once

#include "DiracKernel.h"

#include "GnatBase.h"

class SASFPointSource : public DiracKernel
{
public:
  static InputParameters validParams();

  SASFPointSource(const InputParameters & parameters);

protected:
  virtual void addPoints() override;
  virtual Real computeQpResidual() override;

  // The source location and intensity.
  const Point _source_location;
  const Real _group_source;

  // g
  const unsigned int _group_index;
  // SAAF stabilization parameters.
  const ADMaterialProperty<std::vector<Real>> & _saaf_tau;
}; // class SASFPointSource
