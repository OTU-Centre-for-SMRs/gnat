#pragma once

#include <vector>

#include "Moose.h"
#include "MooseTypes.h"

// The static member function computes a single real spherical harmonics evaluation
// for a given degree l and order m. The class pre-computes all coefficients of
// the complete real spherical harmonics expansion up to a given degree l.
class RealSphericalHarmonics
{
public:
  static Real evaluateCoefficient(unsigned int degree, int order);
  static Real evaluate(unsigned int degree, int order, const Real & mu, const Real & omega);

  RealSphericalHarmonics(unsigned int degree);

  Real evaluatePrecomputed(unsigned int degree, int order, const Real & mu, const Real & omega);

private:
  const unsigned int _degree; // l
  std::vector<Real> _coefficients;
}; // class RealSphericalHarmonics
