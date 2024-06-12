#pragma once

#include "Function.h"

// A function which computes the analytical scalar flux from a point source within annular regions
// (2D flatland) or spherical regions (3D cartesian).
class AxisymmetricPointScalarFlux : public Function
{
public:
  static InputParameters validParams();

  AxisymmetricPointScalarFlux(const InputParameters & parameters);

  virtual Real value(Real t, const Point & p) const override;

protected:
  // The dimensionality of the problem.
  const unsigned int _dims;

  // The point source emission strength.
  const Real _source_strength;
  // The location of the point source.
  const Point _source_location;

  // The radii for the annular axisymmetric regions.
  const std::vector<Real> & _radii;
  // The total macroscopic cross-sections for those regions.
  const std::vector<Real> & _total_xs;
}; // class AxisymmetricPointScalarFlux.
