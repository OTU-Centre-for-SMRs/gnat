#include "RealSphericalHarmonics.h"

#include <cmath>

namespace RealSphericalHarmonics
{
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

/**
 * This hand-coded implementation of the first three degrees of real spherical harmonics
 * is based on the work done by Dr. Park in Moltres:
 * https://github.com/arfc/moltres/blob/devel/src/utils/MoltresUtils.C#L246
 */
Real
evaluate(unsigned int degree, int order, const Real & mu, const Real & eta, const Real & xi)
{
  const Real sqrt_2 = 1.4142135623730951;
  const unsigned int abs_m = static_cast<unsigned int>(std::abs(order));
  const Real C = normalizationConstant(degree, abs_m);

  switch (degree)
  {
    case 0u:
      return 1.0;
    case 1u:
      switch (order)
      {
        case 1:
          return -sqrt_2 * C * eta;
        case 0:
          return C * mu;
        case -1:
          return -sqrt_2 * C * xi;
        default:
          mooseError("Invalid order.");
      }
    case 2:
      switch (order)
      {
        case 2:
          return sqrt_2 * C * 3 * (eta * eta - xi * xi);
        case 1:
          return -sqrt_2 * C * 3 * mu * eta;
        case 0:
          return C * 0.5 * (3 * mu * mu - 1);
        case -1:
          return -sqrt_2 * C * 3 * mu * xi;
        case -2:
          return sqrt_2 * C * 6 * eta * xi;
        default:
          mooseError("Invalid order.");
      }
    case 3:
      switch (order)
      {
        case 3:
          return -sqrt_2 * C * 15 * (4 * eta * eta * eta - 3 * eta * (1 - mu * mu));
        case 2:
          return sqrt_2 * C * 15 * mu * (eta * eta - xi * xi);
        case 1:
          return sqrt_2 * C * 1.5 * (1 - 5 * mu * mu) * eta;
        case 0:
          return C * 0.5 * (5 * mu * mu * mu - 3 * mu);
        case -1:
          return sqrt_2 * C * 1.5 * (1 - 5 * mu * mu) * xi;
        case -2:
          return sqrt_2 * C * 30 * mu * eta * xi;
        case -3:
          return -sqrt_2 * C * 15 * (3 * xi * (1 - mu * mu) - 4 * xi * xi * xi);
        default:
          mooseError("Invalid order.");
      }
    default:
      mooseError("Invalid degree.");
  }
}
} // namespace RealSphericalHarmonics
