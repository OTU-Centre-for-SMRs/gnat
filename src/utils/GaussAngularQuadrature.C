// Implements an angular quadrature set which uses Gauss-Legendre quadrature for
// the polar angular integral and Gauss-Chebyshev quadrature for the azimuthal
// angular integral.
#include "GaussAngularQuadrature.h"

GaussAngularQuadrature::GaussAngularQuadrature(unsigned int n_c,
                                               unsigned int n_l,
                                               MajorAxis axis,
                                               ProblemType type)
  : _n_c(std::move(n_c))
  , _n_l(std::move(n_l))
  , _axis(std::move(axis))
  , _polar_quadrature(std::move(n_l))
  , _azimuthal_quadrature(std::move(n_c))
{
  // Build the quadrature set.
  switch (axis)
  {
    case MajorAxis::X:
    {
      for (unsigned int i = 1; i <= _n_l; ++i)
      {
        for (unsigned int j = 1; j <= _n_c; ++j)
        {
          auto weight = _polar_quadrature.weight(i - 1u)
                        * _azimuthal_quadrature.weight(j - 1u);
          const auto mu = _polar_quadrature.root(i - 1u);
          const auto omega = _azimuthal_quadrature.angularRoot(j - 1u);

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
        }
      }

      break;
    }

    case MajorAxis::Y:
    {
      for (unsigned int i = 1; i <= _n_l; ++i)
      {
        for (unsigned int j = 1; j <= _n_c; ++j)
        {
          auto weight = _polar_quadrature.weight(i - 1u)
                        * _azimuthal_quadrature.weight(j - 1u);
          const auto mu = _polar_quadrature.root(i - 1u);
          const auto omega = _azimuthal_quadrature.angularRoot(j - 1u);

          _quadrature_set_omega.emplace_back(
            RealVectorValue(std::sqrt(1.0 - (mu * mu)) * std::sin(omega), mu,
                            std::sqrt(1.0 - (mu * mu)) * std::cos(omega)));
          _quadrature_set_weight.emplace_back(weight);

          _quadrature_set_omega.emplace_back(
            RealVectorValue(-1.0 * std::sqrt(1.0 - (mu * mu)) * std::sin(omega),
                            -1.0 * mu,
                            -1.0 * std::sqrt(1.0 - (mu * mu)) * std::cos(omega)));
          _quadrature_set_weight.emplace_back(weight);
        }
      }

      break;
    }

    case MajorAxis::Z:
    {
      for (unsigned int i = 1; i <= _n_l; ++i)
      {
        for (unsigned int j = 1; j <= _n_c; ++j)
        {
          auto weight = _polar_quadrature.weight(i - 1u)
                        * _azimuthal_quadrature.weight(j - 1u);
          const auto mu = _polar_quadrature.root(i - 1u);
          const auto omega = _azimuthal_quadrature.angularRoot(j - 1u);

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
        }
      }

      break;
    }
  }

  // Process the quadrature set to remove elements that don't work for certain
  // problem dimensionalities.
  switch (type)
  {
    case ProblemType::Cartesian1D:
      _quadrature_set_omega.clear();
      _quadrature_set_weight.clear();
      for (unsigned int i = 1; i <= _n_l; ++i)
      {
        _quadrature_set_omega.emplace_back(RealVectorValue(_polar_quadrature.root(i - 1u), 0.0, 0.0));
        _quadrature_set_weight.emplace_back(_polar_quadrature.weight(i - 1u));
      }
      break;

    case ProblemType::Cartesian2D:
      for (unsigned int i = 0; i < _quadrature_set_omega.size(); ++i)
      {
        if (_quadrature_set_omega[i](2) < 0.0)
        {
          _quadrature_set_omega.erase(_quadrature_set_omega.begin() + i);
          _quadrature_set_weight.erase(_quadrature_set_weight.begin() + i);
        }
      }
      break;

    default:
      break;
  }
}
