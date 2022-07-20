#pragma once

#include "DiracKernel.h"

#include "GaussAngularQuadrature.h"

#include <map>

// This source is an isotopic neutron source for the current group.
class IsotropicNeutronPointSource : public DiracKernel
{
public:
  static InputParameters validParams();

  IsotropicNeutronPointSource(const InputParameters & parameters);

  virtual void addPoints() override;

protected:
  virtual Real computeQpResidual() override;

  const GaussAngularQuadrature _quadrature_set;
  Real _symmetry_factor;

  const std::vector<Real> _source_intensities; // S_{g, 0, 0} for all points.
  const std::vector<Point> _source_locations;
  std::map<Point, unsigned int> _point_intensity_mapping;
}; // class IsotropicNeutronPointSource
