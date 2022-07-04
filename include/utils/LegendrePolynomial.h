// Implements Legendre polynomials with the intent of computing integrals with
// quadrature sets.
#pragma once

#include <vector>

#include "Moose.h"
#include "MooseTypes.h"

// Implementation inspired by:
// http://rosettacode.org/wiki/Numerical_integration/Gauss-Legendre_Quadrature
class LegendrePolynomial
{
public:
  LegendrePolynomial(unsigned int degree);

  unsigned int degree() const { return _degree; }
  Real root(unsigned int index) const { return _roots[index]; }
  Real weight(unsigned int index) const { return _weights[index]; }

  std::vector<Real> getRoots() const { return _roots; }
  std::vector<Real> getWeights() const { return _weights; }

private:
  // Evaluate the value of a Legendre polynomial at x using the recurrance
  // relationship.
  Real evaluate(Real x);

  // Evaluate the value of a Legendre polynomial's derivative at x using the
  // recurrance relationship.
  Real evaluateDerivative(Real x);

  const unsigned int _degree;

  // Gauss-Legendre weights are easier to compute while finding zeros.
  std::vector<Real> _roots;
  std::vector<Real> _weights;
}; // class LegendrePolynomial
