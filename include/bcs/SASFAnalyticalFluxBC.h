#pragma once

#include "IntegratedBC.h"

class SASFAnalyticalFluxBC : public IntegratedBC
{
public:
  static InputParameters validParams();

  SASFAnalyticalFluxBC(const InputParameters & parameters);

protected:
  virtual Real computeQpResidual() override;

  // Material properties in the near-source region.
  const unsigned int _group_index;
  const std::vector<Real> _sigma_t_g;
  const std::vector<Real> _group_source;

  // The source location.
  const Point _source_location;

  // The number of mesh dimensions.
  const unsigned int _dim;
}; // class SASFAnalyticalFluxBC
