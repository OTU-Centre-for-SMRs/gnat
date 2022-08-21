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
  GaussAngularQuadrature(unsigned int n_c, unsigned int n_l,
                         MajorAxis axis = MajorAxis::X,
                         ProblemType type = ProblemType::Cartesian3D);

  unsigned int totalOrder() const { return _quadrature_set_omega.size(); }
  unsigned int chebyshevOrder() const { return _n_c; }
  unsigned int legendreOrder() const { return _n_l; }
  MajorAxis getAxis() const { return _axis; }
  ProblemType getProblemType() const { return _type; }

  RealVectorValue direction(unsigned int n) const { return _quadrature_set_omega[n]; }
  Real weight(unsigned int n) const { return _quadrature_set_weight[n]; }

  std::vector<RealVectorValue> & getDirections() { return _quadrature_set_omega; }
  std::vector<Real> & getWeights() { return _quadrature_set_weight; }

  LegendrePolynomial getPolarLegendre() const { return _polar_quadrature; }
  ChebyshevPolynomial getAzimuthalChebyshev() const { return _azimuthal_quadrature; }

private:
  // Generate a weight-ordinate pair for 1-3 dimensional problems.
  void generateWeightOrdiantePair(const unsigned int i, const unsigned int j);

  // Number of Chebyshev quadrature points and number of Legendre quadrature
  // points, respectively.
  const unsigned int _n_c;
  const unsigned int _n_l;
  const MajorAxis _axis;
  const ProblemType _type;

  LegendrePolynomial _polar_quadrature;
  ChebyshevPolynomial _azimuthal_quadrature;

  std::vector<RealVectorValue> _quadrature_set_omega;
  std::vector<Real> _quadrature_set_weight;
}; // class GaussAngularQuadrature
