// Implements an angular quadrature set which uses Gauss-Legendre quadrature for
// the polar angular integral and Gauss-Chebyshev quadrature for the azimuthal
// angular integral.
#include "GaussAngularQuadrature.h"

GaussAngularQuadrature::GaussAngularQuadrature(unsigned int n_c,
                                               unsigned int n_l,
                                               MajorAxis axis,
                                               ProblemType type)
  : AngularQuadrature(std::move(axis), std::move(type)),
    _n_c(std::move(n_c)),
    _n_l(std::move(n_l)),
    _polar_quadrature(std::move(n_l)),
    _azimuthal_quadrature(std::move(n_c))
{
  switch (type)
  {
    case ProblemType::Cartesian1D:
      for (unsigned int i = 1u; i <= _n_l; ++i)
        generateWeightOrdiantePair(i, 0u);

      break;

    case ProblemType::Cartesian2D:
      for (unsigned int i = 1u; i <= _n_l; ++i)
      {
        for (unsigned int j = 1u; j <= _n_c; ++j)
          generateWeightOrdiantePair(i, j);
      }

      break;

    case ProblemType::Cartesian3D:
      for (unsigned int i = 1u; i <= _n_l; ++i)
      {
        for (unsigned int j = 1u; j <= _n_c; ++j)
          generateWeightOrdiantePair(i, j);
      }

      break;
  }
}

unsigned int
GaussAngularQuadrature::totalOrder() const
{
  return _quadrature_set_omega.size();
}
const RealVectorValue &
GaussAngularQuadrature::direction(unsigned int n) const
{
  return _quadrature_set_omega[n];
}
const Real &
GaussAngularQuadrature::weight(unsigned int n) const
{
  return _quadrature_set_weight[n];
}
const std::vector<RealVectorValue> &
GaussAngularQuadrature::getDirections() const
{
  return _quadrature_set_omega;
}
const std::vector<Real> &
GaussAngularQuadrature::getWeights() const
{
  return _quadrature_set_weight;
}

const Real &
GaussAngularQuadrature::getPolarRoot(unsigned int n) const
{
  return _polar_quadrature.root(n / _n_l);
}
const Real &
GaussAngularQuadrature::getAzimuthalAngularRoot(unsigned int n) const
{
  return _azimuthal_quadrature.angularRoot(n % _n_c);
}

void
GaussAngularQuadrature::generateWeightOrdiantePair(unsigned int i, unsigned int j)
{
  Real weight = 0.0;
  Real mu = 0.0;
  Real omega = 0.0;

  switch (_type)
  {
    case ProblemType::Cartesian1D:
      _quadrature_set_omega.emplace_back(RealVectorValue(_polar_quadrature.root(i - 1u), 0.0, 0.0));
      _quadrature_set_weight.emplace_back(_polar_quadrature.weight(i - 1u));
      break;

    case ProblemType::Cartesian2D:
      weight = _polar_quadrature.weight(i - 1u) * _azimuthal_quadrature.weight(j - 1u);
      mu = _polar_quadrature.root(i - 1u);
      omega = _azimuthal_quadrature.angularRoot(j - 1u);
      _quadrature_set_omega.emplace_back(
          RealVectorValue(mu, std::sqrt(1.0 - (mu * mu)) * std::cos(omega), 0.0));
      _quadrature_set_weight.emplace_back(weight);

      break;

    case ProblemType::Cartesian3D:
      // Allow the use of a different major axis for the quadrature set.
      switch (_axis)
      {
        case MajorAxis::X:
          weight = _polar_quadrature.weight(i - 1u) * _azimuthal_quadrature.weight(j - 1u);
          mu = _polar_quadrature.root(i - 1u);
          omega = _azimuthal_quadrature.angularRoot(j - 1u);
          _quadrature_set_omega.emplace_back(
              RealVectorValue(mu,
                              std::sqrt(1.0 - (mu * mu)) * std::cos(omega),
                              std::sqrt(1.0 - (mu * mu)) * std::sin(omega)));
          _quadrature_set_weight.emplace_back(weight);

          _quadrature_set_omega.emplace_back(
              RealVectorValue(-1.0 * mu,
                              -1.0 * std::sqrt(1.0 - (mu * mu)) * std::cos(omega),
                              -1.0 * std::sqrt(1.0 - (mu * mu)) * std::sin(omega)));
          _quadrature_set_weight.emplace_back(weight);

          break;

        case MajorAxis::Y:
          weight = _polar_quadrature.weight(i - 1u) * _azimuthal_quadrature.weight(j - 1u);
          mu = _polar_quadrature.root(i - 1u);
          omega = _azimuthal_quadrature.angularRoot(j - 1u);
          _quadrature_set_omega.emplace_back(
              RealVectorValue(std::sqrt(1.0 - (mu * mu)) * std::sin(omega),
                              mu,
                              std::sqrt(1.0 - (mu * mu)) * std::cos(omega)));
          _quadrature_set_weight.emplace_back(weight);

          _quadrature_set_omega.emplace_back(
              RealVectorValue(-1.0 * std::sqrt(1.0 - (mu * mu)) * std::sin(omega),
                              -1.0 * mu,
                              -1.0 * std::sqrt(1.0 - (mu * mu)) * std::cos(omega)));
          _quadrature_set_weight.emplace_back(weight);

          break;

        case MajorAxis::Z:
          weight = _polar_quadrature.weight(i - 1u) * _azimuthal_quadrature.weight(j - 1u);
          mu = _polar_quadrature.root(i - 1u);
          omega = _azimuthal_quadrature.angularRoot(j - 1u);
          _quadrature_set_omega.emplace_back(
              RealVectorValue(std::sqrt(1.0 - (mu * mu)) * std::cos(omega),
                              std::sqrt(1.0 - (mu * mu)) * std::sin(omega),
                              mu));
          _quadrature_set_weight.emplace_back(weight);

          _quadrature_set_omega.emplace_back(
              RealVectorValue(-1.0 * std::sqrt(1.0 - (mu * mu)) * std::cos(omega),
                              -1.0 * std::sqrt(1.0 - (mu * mu)) * std::sin(omega),
                              -1.0 * mu));
          _quadrature_set_weight.emplace_back(weight);

          break;
      }

      break;
  }
}
