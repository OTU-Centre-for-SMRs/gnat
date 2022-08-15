#include "ADDFEMScattering.h"

#include "RealSphericalHarmonics.h"

registerMooseObject("GnatApp", ADDFEMScattering);

InputParameters
ADDFEMScattering::validParams()
{
  auto params = ADSNBaseKernel::validParams();
  params.addClassDescription("Computes the scattering term for the "
                             "current group of the discrete ordinates neutron "
                             "transport equation. The weak form is given by "
                             "$-(\\psi_{j}, \\sum_{g' = 1}^{G}"
                             "\\Sigma_{s,\\, g'\\rightarrow g}"
                             "\\sum_{l = 0}^{L}\\frac{2l + 1}{4\\pi} "
                             "f_{g'\\rightarrow g,\\, l}"
                             "\\sum_{m = -l}^{l}Y_{l,m}(\\hat{\\Omega}_{n})"
                             "\\Phi_{g',l,m})$. The group flux "
                             "moments are computed by this kernel. This kernel "
                             "should not be exposed to the user, instead being "
                             "enabled through a transport action. This kernel "
                             "is provided for debugging purposes only: "
                             "full-matrix scattering evaluation without source "
                             "iteration is quite slow and should not be used "
                             "for production calculations.");
  params.addRequiredCoupledVar("group_flux_ordinates",
                               "The angular flux ordinates for all groups.");
  params.addRequiredRangeCheckedParam<unsigned int>("ordinate_index",
                                                    "ordinate_index >= 0",
                                                    "The discrete ordinate index "
                                                    "of the current angular "
                                                    "flux.");
  params.addRequiredRangeCheckedParam<unsigned int>("group_index",
                                                    "group_index >= 0",
                                                    "The energy group index "
                                                    "of the current angular "
                                                    "flux.");
  params.addRequiredRangeCheckedParam<unsigned int>("num_groups",
                                                    "num_groups >= 1",
                                                    "The number of spectral "
                                                    "energy groups.");
  params.addRequiredRangeCheckedParam<unsigned int>("max_anisotropy",
                                                    "max_anisotropy >= 0",
                                                    "The maximum degree of "
                                                    "anisotropy to evaluate.");

  return params;
}

ADDFEMScattering::ADDFEMScattering(const InputParameters & parameters)
  : ADSNBaseKernel(parameters)
  , _sigma_s_g_prime_g_l(getADMaterialProperty<std::vector<Real>>("scattering_matrix"))
  , _anisotropy(getMaterialProperty<unsigned int>("medium_anisotropy"))
  , _ordinate_index(getParam<unsigned int>("ordinate_index"))
  , _group_index(getParam<unsigned int>("group_index"))
  , _num_groups(getParam<unsigned int>("num_groups"))
  , _max_anisotropy(getParam<unsigned int>("max_anisotropy"))
{
  if (_group_index >= _num_groups)
    mooseError("The group index exceeds the number of energy groups.");

  if (_ordinate_index >= _quadrature_set.totalOrder())
    mooseError("The ordinates index exceeds the number of quadrature points.");

  const unsigned int num_coupled = coupledComponents("group_flux_ordinates");
  if (num_coupled != _quadrature_set.totalOrder() * _num_groups)
    mooseError("Mismatch between the angular flux ordinates and quadrature set.");

  _group_flux_ordinates.reserve(num_coupled);
  for (unsigned int i = 0; i < num_coupled; ++i)
    _group_flux_ordinates.emplace_back(&coupledValue("group_flux_ordinates", i));
}

ADReal
ADDFEMScattering::computeFluxMoment(unsigned int g_prime, unsigned int l,
                                    int m)
{
  ADReal moment = 0.0;
  Real mu = 0.0;
  Real omega = 0.0;
  for (unsigned int i = 0; i < _quadrature_set.totalOrder(); ++i)
  {
    cartesianToSpherical(_quadrature_set.direction(i),
                         mu, omega);

    const unsigned int base_n = g_prime * (_group_flux_ordinates.size()
                                           / _num_groups);
    moment += RealSphericalHarmonics::evaluate(l, m, mu, omega)
              * std::max(MetaPhysicL::raw_value((* _group_flux_ordinates[base_n + i])[_qp]), 0.0)
              * _quadrature_set.weight(i);
  }

  return moment;
}

// Compute the full scattering term for both in-group and group-to-group
// scattering.
ADReal
ADDFEMScattering::computeQpResidual()
{
  // Quit early if no Legendre cross-section moments are provided.
  if (_sigma_s_g_prime_g_l[_qp].size() == 0u)
    return 0.0;

  ADReal moment_l = 0.0;
  ADReal res = 0.0;
  Real omega = 0.0;
  Real mu = 0.0;
  unsigned int scattering_index = 0u;
  for (unsigned int g_prime = 0; g_prime < _num_groups; ++g_prime)
  {
    scattering_index = g_prime * _num_groups * _anisotropy[_qp]
                       + _group_index * _anisotropy[_qp];

    // The maximum degree of anisotropy we can handle.
    const unsigned int max_anisotropy = std::min(_anisotropy[_qp],
                                                 max_anisotropy);
    for (unsigned int l = 0; l <= max_anisotropy; ++l)
    {
      // Handle different levels of dimensionality.
      switch (_quadrature_set.getProblemType())
      {
        // Legendre moments in 1D, looping over m is unecessary.
        case ProblemType::Cartesian1D:
          cartesianToSpherical(_quadrature_set.direction(_ordinate_index),
                               mu, omega);
          moment_l += computeFluxMoment(g_prime, l, 0)
                      * RealSphericalHarmonics::evaluate(l, 0, mu, omega);
          break;

        // Need moments with m >= 0 for 2D.
        case ProblemType::Cartesian2D:
          for (int m = 0; m <= static_cast<int>(l); ++m)
          {
            cartesianToSpherical(_quadrature_set.direction(_ordinate_index),
                                 mu, omega);
            moment_l += computeFluxMoment(g_prime, l, m)
                        * RealSphericalHarmonics::evaluate(l, m, mu, omega);
          }
          break;

        // Need all moments in 3D.
        case ProblemType::Cartesian3D:
          for (int m = -1 * static_cast<int>(l); m <= static_cast<int>(l); ++m)
          {
            cartesianToSpherical(_quadrature_set.direction(_ordinate_index),
                                 mu, omega);
            moment_l += computeFluxMoment(g_prime, l, m)
                        * RealSphericalHarmonics::evaluate(l, m, mu, omega);
          }
          break;

        default: // Defaults to doing nothing for now.
          break;
      }

      res += (2.0 * static_cast<Real>(l) + 1.0) / (4.0 * M_PI)
             * _sigma_s_g_prime_g_l[_qp][scattering_index] * moment_l
             * _symmetry_factor;

      scattering_index++;
      moment_l = 0.0;
    }
  }

  return -1.0 * _test[_i][_qp] * res;
}
