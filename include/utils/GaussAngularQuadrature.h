// Implements an angular quadrature set which uses Gauss-Legendre quadrature for
// the polar angular dependance and Gauss-Chebyshev quadrature for the azimuthal
// angular dependance.
#pragma once

#include <vector>

#include "MooseTypes.h"
#include "libmesh/vector_value.h"

#include "LegendrePolynomial.h"
#include "ChebyshevPolynomial.h"

#include "AngularQuadrature.h"

class GaussAngularQuadrature : public AngularQuadrature
{
public:
  GaussAngularQuadrature(unsigned int n_c,
                         unsigned int n_l,
                         MajorAxis axis = MajorAxis::X,
                         ProblemType type = ProblemType::Cartesian3D);

  unsigned int totalOrder() const override;
  const RealVectorValue & direction(unsigned int n) const override;
  const Real & weight(unsigned int n) const override;
  const std::vector<RealVectorValue> & getDirections() const override;
  const std::vector<Real> & getWeights() const override;

  const Real & getPolarRoot(unsigned int n) const override;
  const Real & getAzimuthalAngularRoot(unsigned int n) const override;

  unsigned int legendreOrder() const { return _n_l; }
  LegendrePolynomial getPolarLegendre() const { return _polar_quadrature; }

  unsigned int chebyshevOrder() const { return _n_c; }
  ChebyshevPolynomial getAzimuthalChebyshev() const { return _azimuthal_quadrature; }

private:
  // Generate a weight-ordinate pair for 1-3 dimensional problems.
  void generateWeightOrdiantePair(unsigned int i, unsigned int j);

  // Number of Chebyshev quadrature points and number of Legendre quadrature
  // points, respectively.
  const unsigned int _n_c;
  const unsigned int _n_l;

  LegendrePolynomial _polar_quadrature;
  ChebyshevPolynomial _azimuthal_quadrature;

  std::vector<RealVectorValue> _quadrature_set_omega;
  std::vector<Real> _quadrature_set_weight;
}; // class GaussAngularQuadrature
