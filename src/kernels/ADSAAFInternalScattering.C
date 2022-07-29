#include "ADSAAFInternalScattering.h"

registerMooseObject("GnatApp", ADSAAFInternalScattering);

InputParameters
ADSAAFInternalScattering::validParams()
{
  auto params = ADSAAFBaseKernel::validParams();
  params.addClassDescription("Computes the within-group scattering source for the "
                             "current group of the SAAF discrete ordinates neutron "
                             "transport equation. The weak form is given by "
                             "$-(\\psi_{j} + \\tau_{g}\\vec{\\nabla}\\psi_{j}"
                             "\\cdot\\hat{\\Omega}, \\Sigma_{s,\\, g\\rightarrow g}"
                             "\\sum_{l = 0}^{L}\\frac{2l + 1}{4\\pi} "
                             "f_{g\\rightarrow g,\\, l}"
                             "\\sum_{m = -l}^{l}Y_{l,m}(\\hat{\\Omega}_{n})"
                             "\\Phi_{g,l,m})$. The flux moments required for "
                             "this kernel are expected to be supplied as "
                             "Auxvariables. This kernel should not be exposed "
                             "to the user, instead being enabled through a "
                             "transport action.");
  params.addRequiredCoupledVar("within_group_flux_moments",
                               "The flux moments of the current group.");
  params.addRequiredRangeCheckedParam<unsigned int>("num_groups",
                                                    "num_groups >= 1",
                                                    "The number of spectral "
                                                    "energy groups.");

  return params;
}

ADSAAFInternalScattering::ADSAAFInternalScattering(const InputParameters & parameters)
  : ADSAAFBaseKernel(parameters)
  , _sigma_s_g_prime_g_l(getADMaterialProperty<std::vector<Real>>("scattering_matrix"))
  , _anisotropy(getMaterialProperty<unsigned int>("medium_anisotropy"))
  , _num_groups(getParam<unsigned int>("num_groups"))
{
  if (_group_index >= _num_groups)
    mooseError("The group index exceeds the number of energy groups.");

  if (_ordinate_index >= _quadrature_set.totalOrder())
    mooseError("The ordinates index exceeds the number of quadrature points.");

  const unsigned int num_coupled = coupledComponents("within_group_flux_moments");

  switch (_quadrature_set.getProblemType())
  {
    case ProblemType::Cartesian1D:
      _provided_moment_degree = num_coupled - 1u;
      break;

    case ProblemType::Cartesian2D:
      // Doing it this way to minimize divisions. Trying to keep it fixed-point
      // for as long as possible.
      _provided_moment_degree = static_cast<unsigned int>(-3
                                  + std::sqrt(9 - 8 * (1
                                    - static_cast<int>(num_coupled))));
      _provided_moment_degree /= 2u;
      break;

    case ProblemType::Cartesian3D:
      _provided_moment_degree = std::sqrt(num_coupled) - 1u;
      break;

    default:
      _provided_moment_degree = 0u;

      _symmetry_factor = 1.0;
      break;
  }

  _within_group_flux_moments.reserve(num_coupled);
  for (unsigned int i = 0; i < num_coupled; ++i)
    _within_group_flux_moments.emplace_back(&coupledValue("within_group_flux_moments", i));
}

ADReal
ADSAAFInternalScattering::computeQpResidual()
{
  // Quit early if no Legendre cross-section moments are provided.
  if (_sigma_s_g_prime_g_l[_qp].size() == 0u)
    return 0.0;

  ADReal res, moment_l = 0.0;
  Real omega, mu = 0.0;

  const unsigned int scattering_index = _group_index * _num_groups
                                        * _anisotropy[_qp]
                                        + _group_index * _anisotropy[_qp];
  unsigned int moment_index = 0u;
  // The maximum degree of anisotropy we can handle.
  const unsigned int max_anisotropy = std::min(_anisotropy[_qp],
                                               _provided_moment_degree);
  for (unsigned int l = 0; l <= max_anisotropy; ++l)
  {
    // Handle different levels of dimensionality.
    switch (_quadrature_set.getProblemType())
    {
      // Legendre moments in 1D, looping over m is unecessary.
      case ProblemType::Cartesian1D:
        cartesianToSpherical(_quadrature_set.direction(_ordinate_index),
                             mu, omega);
        moment_l += (*_within_group_flux_moments[moment_index])[_qp]
                    * RealSphericalHarmonics::evaluate(l, 0, mu, omega);
        moment_index++;
        break;

      // Need moments with m >= 0 for 2D.
      case ProblemType::Cartesian2D:
        for (int m = 0; m <= static_cast<int>(l); ++m)
        {
          cartesianToSpherical(_quadrature_set.direction(_ordinate_index),
                               mu, omega);
          moment_l += (*_within_group_flux_moments[moment_index])[_qp]
                      * RealSphericalHarmonics::evaluate(l, m, mu, omega);
          moment_index++;
        }
        break;

      // Need all moments in 3D.
      case ProblemType::Cartesian3D:
        for (int m = -1 * static_cast<int>(l); m <= static_cast<int>(l); ++m)
        {
          cartesianToSpherical(_quadrature_set.direction(_ordinate_index),
                               mu, omega);
          moment_l += (*_within_group_flux_moments[moment_index])[_qp]
                      * RealSphericalHarmonics::evaluate(l, m, mu, omega);
          moment_index++;
        }
        break;

      default: // Defaults to doing nothing for now.
        break;
    }

    res += (2.0 * static_cast<Real>(l) + 1.0) / (4.0 * M_PI)
           * _sigma_s_g_prime_g_l[_qp][scattering_index + l] * moment_l
           * _symmetry_factor;
    moment_l = 0.0;
  }

  return -1.0 * computeQPTests() * res;
}
