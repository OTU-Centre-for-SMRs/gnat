// Implements Chebyshev polynomials with the intent of computing integrals with
// quadrature sets.
#pragma once

#include <vector>

#include "Moose.h"
#include "MooseTypes.h"

class ChebyshevPolynomial
{
public:
  ChebyshevPolynomial(unsigned int degree);

  unsigned int degree() const { return _degree; }
  const Real & root(unsigned int index) const { return _roots[index]; }
  const Real & angularRoot(unsigned int index) const { return _angular_roots[index]; }
  const Real & weight(unsigned int index) const { return _weights[index]; }

  const std::vector<Real> & getRoots() const { return _roots; }
  const std::vector<Real> & getAngularRoots() const { return _angular_roots; }
  const std::vector<Real> & getWeights() const { return _weights; }

private:
  const unsigned int _degree;

  // Storing weights to mirror LegendrePolynomial.
  // _roots stores the roots on -1 \leq y \leq 1
  // _angular_roots stores the roots on 0 \leq \omega \leq 2\pi
  std::vector<Real> _roots;
  std::vector<Real> _angular_roots;
  std::vector<Real> _weights;
}; // class ChebyshevPolynomial
