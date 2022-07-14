// Implements an angular quadrature set which uses Gauss-Legendre quadrature for
// the polar angular dependance and Gauss-Chebyshev quadrature for the azimuthal
// angular dependance.
#pragma once

#include <vector>

#include "Moose.h"
#include "MooseTypes.h"
#include "libmesh/vector_value.h"

#include "LegendrePolynomial.h"
#include "ChebyshevPolynomial.h"

#include "GnatBase.h"

class GaussAngularQuadrature
{
public:
  GaussAngularQuadrature(unsigned int n_c, unsigned int n_l, MajorAxis axis = MajorAxis::X);

  unsigned int totalOrder() const { return 2 * _n_c * _n_l; }
  unsigned int chebyshevOrder() const { return _n_c; }
  unsigned int legendreOrder() const { return _n_l; }
  MajorAxis getAxis() const { return _axis; }

  RealVectorValue getQPDirection(unsigned int n) const { return _quadrature_set_omega[n]; }
  Real getQPWeight(unsigned int n) const { return _quadrature_set_weight[n]; }

  std::vector<RealVectorValue> & getDirections() { return _quadrature_set_omega; }
  std::vector<Real> & getWeights() { return _quadrature_set_weight; }

  LegendrePolynomial getPolarLegendre() const { return _polar_quadrature; }
  ChebyshevPolynomial getAzimuthalChebyshev() const { return _azimuthal_quadrature; }

private:
  // Number of Chebyshev quadrature points and number of Legendre quadrature
  // points, respectively.
  const unsigned int _n_c;
  const unsigned int _n_l;
  const MajorAxis _axis;

  LegendrePolynomial _polar_quadrature;
  ChebyshevPolynomial _azimuthal_quadrature;

  std::vector<RealVectorValue> _quadrature_set_omega;
  std::vector<Real> _quadrature_set_weight;
}; // class GaussAngularQuadrature
