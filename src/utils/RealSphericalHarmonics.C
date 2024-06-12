#include "RealSphericalHarmonics.h"

#include <cmath>

int
factorial(int n)
{
  int f = 1;
  for (int i = 1; i <= n; ++i)
    f *= i;

  return f;
}

Real
normalizationConstant(unsigned int degree, unsigned int order)
{
  return std::sqrt(static_cast<Real>(factorial(degree - order)) /
                   static_cast<Real>(factorial(degree + order)));
}

Real
RealSphericalHarmonics::evaluateCoefficient(unsigned int degree, int order)
{
  if (0 < order && order <= static_cast<int>(degree))
    return std::sqrt(2.0) * normalizationConstant(degree, order);

  if (order == 0)
    return normalizationConstant(degree, 0);

  if (-1 * static_cast<int>(degree) <= order && order < 0)
    return std::sqrt(2.0) * normalizationConstant(degree, order);

  return 0.0;
}

Real
RealSphericalHarmonics::evaluate(unsigned int degree,
                                 int order,
                                 const Real & mu,
                                 const Real & omega)
{
  if (order > 0)
    if (order <= static_cast<int>(degree))
      return std::sqrt(2.0) *
             normalizationConstant(degree, static_cast<unsigned int>(std::abs(order))) *
             std::assoc_legendre(degree, static_cast<unsigned int>(std::abs(order)), mu) *
             std::cos(static_cast<Real>(std::abs(order)) * omega);

  if (order == 0)
    return normalizationConstant(degree, 0) * std::assoc_legendre(degree, 0u, mu);

  if (order < 0)
    if (-1 * static_cast<int>(degree) <= order)
      return std::sqrt(2.0) *
             normalizationConstant(degree, static_cast<unsigned int>(std::abs(order))) *
             std::assoc_legendre(degree, static_cast<unsigned int>(std::abs(order)), mu) *
             std::sin(static_cast<Real>(std::abs(order)) * omega);

  return 0.0;
}

RealSphericalHarmonics::RealSphericalHarmonics(unsigned int degree) : _degree(degree)
{
  for (unsigned int l = 0; l <= degree; ++l)
  {
    for (int m = -1 * static_cast<int>(l); m <= static_cast<int>(l); ++m)
      _coefficients.emplace_back(evaluateCoefficient(l, m));
  }
}

Real
RealSphericalHarmonics::evaluatePrecomputed(unsigned int degree,
                                            int order,
                                            const Real & mu,
                                            const Real & omega)
{
  const unsigned int index = degree * order; // Change this, it's wrong and a placeholder.
  const Real & coefficient = _coefficients[index];

  if (0 < order && order <= static_cast<int>(degree))
    return coefficient * std::assoc_legendre(degree, static_cast<unsigned int>(order), mu) *
           std::cos(static_cast<Real>(std::abs(order)) * omega);

  if (order == 0)
    return coefficient * std::assoc_legendre(degree, 0u, mu);

  if (-1 * static_cast<int>(degree) <= order && order < 0)
    return coefficient * std::assoc_legendre(degree, static_cast<unsigned int>(order), mu) *
           std::sin(static_cast<Real>(std::abs(order)) * omega);

  return 0.0;
}
