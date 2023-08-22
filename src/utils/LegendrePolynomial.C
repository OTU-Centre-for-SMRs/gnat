// Implements Legendre polynomials with the intent of computing integrals with
// quadrature sets.
#include "LegendrePolynomial.h"

#include <cmath>

//------------------------------------------------------------------------------
// Legendre polynomials here.
//------------------------------------------------------------------------------
LegendrePolynomial::LegendrePolynomial(unsigned int degree) : _degree(std::move(degree))
{
  // Precompute and store the roots of the polynomial using Newton's method.
  _roots.resize(degree, 0.0);
  _weights.resize(degree, 0.0);

  Real dr = 0.0;
  Real x = 0.0;
  Real v = 0.0;
  Real d = 0.0;
  for (unsigned int i = 1; i <= _degree; ++i)
  {
    dr = 1;

    x = std::cos(libMesh::pi * (static_cast<Real>(i - 0.25)) /
                 (static_cast<Real>(_degree) + 0.5)); // Magical initial guess.
    v = evaluate(x);
    d = evaluateDerivative(x);

    do
    {
      dr = v / d;
      x -= dr; // x_{i+1} = x_{i} - f(x_{i}) / f'(x_{i})
      v = evaluate(x);
      d = evaluateDerivative(x);
    } while (std::abs(dr) > 2e-16);

    _roots[i - 1u] = x;
    _weights[i - 1u] = 2.0 / ((1.0 - x * x) * d * d);
  }
}

// Evaluate the value of a Legendre polynomial at x using the recurrance
// relationship.
Real
LegendrePolynomial::evaluate(Real x)
{
  if (_degree == 0u)
    return 1.0;
  if (_degree == 1u)
    return x;

  Real v = 1.0;
  Real v_sub_1 = x;
  Real v_sub_2 = 1.0;

  for (unsigned int i = 2; i <= _degree; ++i)
  {
    v = ((2.0 * static_cast<Real>(i) - 1.0) * x * v_sub_1 - (i - 1.0) * v_sub_2) /
        static_cast<Real>(i);

    v_sub_2 = v_sub_1;
    v_sub_1 = v;
  }

  return v;
}

// Evaluate the value of a Legendre polynomial's derivative at x using the
// recurrance relationship.
Real
LegendrePolynomial::evaluateDerivative(Real x)
{
  if (_degree == 0u)
    return 0.0;
  if (_degree == 1u)
    return 1.0;

  Real v = 1.0;
  Real d = 0.0;
  Real v_sub_1 = x;
  Real v_sub_2 = 1.0;
  Real f = 1.0 / (x * x - 1.0);

  for (unsigned int i = 2; i <= _degree; ++i)
  {
    v = ((2.0 * static_cast<Real>(i) - 1.0) * x * v_sub_1 -
         (static_cast<Real>(i) - 1.0) * v_sub_2) /
        static_cast<Real>(i);
    d = static_cast<Real>(i) * f * (x * v - v_sub_1);

    v_sub_2 = v_sub_1;
    v_sub_1 = v;
  }

  return d;
}
