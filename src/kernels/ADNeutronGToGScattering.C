#include "ADNeutronGToGScattering.h"

#include "RealSphericalHarmonics.h"

registerMooseObject("GnatApp", ADNeutronGToGScattering);

InputParameters
ADNeutronGToGScattering::validParams()
{
  auto params = ADKernel::validParams();
  params.addClassDescription("Computes the group-to-group scattering term for the "
                             "current group of the discrete ordinates neutron "
                             "transport equation. The weak form is given by "
                             "$-(\\psi_{j}, \\sum_{g' = 1,\\, g'\\neq g}^{G}"
                             "\\Sigma_{s,\\, g'\\rightarrow g}"
                             "\\sum_{l = 0}^{L}\\frac{2l + 1}{4\\pi} "
                             "f_{g'\\rightarrow g,\\, l}"
                             "\\sum_{m = -l}^{l}Y_{l,m}(\\hat{\\Omega}_{n})"
                             "\\Phi_{g',l,m})$. The group flux "
                             "moments required for this kernel are expected to "
                             "be supplied as Auxvariables. "
                             "This kernel should not be exposed to the user, "
                             "instead being enabled through a transport action.");
  params.addRequiredCoupledVar("group_flux_moments",
                               "The flux moments for all groups.");
  params.addRequiredParam<MooseEnum>("dimensionality",
                                     MooseEnum("1D_cartesian 2D_cartesian 3D_cartesian"),
                                     "Dimensionality and the coordinate system of the "
                                     "problem.");
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

  return params;
}

ADNeutronGToGScattering::ADNeutronGToGScattering(const InputParameters & parameters)
  : ADKernel(parameters)
  , _type(getParam<MooseEnum>("dimensionality").getEnum<ProblemType>())
  , _directions(getMaterialProperty<std::vector<RealVectorValue>>("directions"))
  , _axis(getMaterialProperty<MajorAxis>("quadrature_axis_alignment"))
  , _sigma_s_g_prime_g_l(getADMaterialProperty<std::vector<Real>>("scattering_matrix"))
  , _anisotropy(getMaterialProperty<unsigned int>("medium_anisotropy"))
  , _ordinate_index(getParam<unsigned int>("ordinate_index"))
  , _group_index(getParam<unsigned int>("group_index"))
  , _num_groups(getParam<unsigned int>("num_groups"))
  , _symmetry_factor(1.0)
{
  if (_group_index >= _num_groups)
    mooseError("The group index exceeds the number of energy groups.");

  const unsigned int num_coupled = coupledComponents("group_flux_moments");
  const unsigned int num_group_moments = num_coupled / _num_groups;

  switch (_type)
  {
    case ProblemType::Cartesian1D:
      _provided_moment_degree = num_group_moments - 1u;

      _symmetry_factor = 2.0 * M_PI;
      break;

    case ProblemType::Cartesian2D:
      // Doing it this way to minimize divisions. Trying to keep it fixed-point
      // for as long as possible.
      _provided_moment_degree = static_cast<unsigned int>(-3 + std::sqrt(9 - 8 * (1 - static_cast<int>(num_group_moments))));
      _provided_moment_degree /= 2u;

      _symmetry_factor = 2.0;
      break;

    case ProblemType::Cartesian3D:
      _provided_moment_degree = std::sqrt(num_group_moments) - 1u;

      _symmetry_factor = 1.0;
      break;

    default:
      _provided_moment_degree = 0u;

      _symmetry_factor = 1.0;
      break;
  }

  _group_flux_moments.reserve(num_coupled);
  for (unsigned int i = 0; i < num_coupled; ++i)
    _group_flux_moments.emplace_back(&coupledValue("group_flux_moments", i));
}

void
ADNeutronGToGScattering::cartesianToSpherical(const RealVectorValue & ordinate,
                                              Real & mu, Real & omega)
{
  switch (_axis[_qp])
  {
    case MajorAxis::X:
      mu = ordinate(0);
      omega = std::acos(ordinate(1) / std::sqrt(1.0 - (mu * mu)));

      break;

    case MajorAxis::Y:
      mu = ordinate(1);
      omega = std::acos(ordinate(2) / std::sqrt(1.0 - (mu * mu)));

      break;

    case MajorAxis::Z:
      mu = ordinate(2);
      omega = std::acos(ordinate(0) / std::sqrt(1.0 - (mu * mu)));

      break;
  }
}

ADReal
ADNeutronGToGScattering::computeQpResidual()
{
  if (_ordinate_index >= _directions[_qp].size())
    mooseError("The ordinates index exceeds the number of quadrature points.");

  // Quit early if no Legendre cross-section moments are provided.
  if (_sigma_s_g_prime_g_l[_qp].size() == 0u)
    return 0.0;

  ADReal res, moment_l = 0.0;
  Real omega, mu = 0.0;

  unsigned int num_moments;
  switch (_type)
  {
    case ProblemType::Cartesian1D:
      num_moments = _provided_moment_degree + 1u;
      break;

    case ProblemType::Cartesian2D:
      num_moments = (_provided_moment_degree + 1u)
                    * (_provided_moment_degree + 2u) / 2u;
      break;

    case ProblemType::Cartesian3D:
      num_moments = (_provided_moment_degree + 1u)
                    * (_provided_moment_degree + 1u);
      break;

    default:
      num_moments = (_provided_moment_degree + 1u)
                    * (_provided_moment_degree + 1u);
      break;
  }

  unsigned int scattering_index = 0u;
  unsigned int moment_index = 0u;
  for (unsigned int g_prime = 0; g_prime < _num_groups; ++g_prime)
  {
    moment_index = g_prime * num_moments;
    scattering_index = g_prime * _num_groups * _anisotropy[_qp]
                       + _group_index * _anisotropy[_qp];

    if (g_prime == _group_index)
      continue;

    // The maximum degree of anisotropy we can handle.
    const unsigned int max_anisotropy = std::min(_anisotropy[_qp],
                                                 _provided_moment_degree);
    for (unsigned int l = 0; l <= max_anisotropy; ++l)
    {
      // Handle different levels of dimensionality.
      switch (_type)
      {
        // Legendre moments in 1D, looping over m is unecessary.
        case ProblemType::Cartesian1D:
          cartesianToSpherical(MetaPhysicL::raw_value(_directions[_qp][_ordinate_index]),
                               mu, omega);
          moment_l += (*_group_flux_moments[moment_index])[_qp]
                      * RealSphericalHarmonics::evaluate(l, 0, mu, omega);
          moment_index++;
          break;

        // Need moments with m >= 0 for 2D.
        case ProblemType::Cartesian2D:
          for (int m = 0; m <= static_cast<int>(l); ++m)
          {
            cartesianToSpherical(MetaPhysicL::raw_value(_directions[_qp][_ordinate_index]),
                                 mu, omega);
            moment_l += (*_group_flux_moments[moment_index])[_qp]
                        * RealSphericalHarmonics::evaluate(l, m, mu, omega);
            moment_index++;
          }
          break;

        // Need all moments in 3D.
        case ProblemType::Cartesian3D:
          for (int m = -1 * static_cast<int>(l); m <= static_cast<int>(l); ++m)
          {
            cartesianToSpherical(MetaPhysicL::raw_value(_directions[_qp][_ordinate_index]),
                                 mu, omega);
            moment_l += (*_group_flux_moments[moment_index])[_qp]
                        * RealSphericalHarmonics::evaluate(l, m, mu, omega);
            moment_index++;
          }
          break;

        default: // Defaults to doing nothing for now.
          break;
      }

      res += (2.0 * static_cast<Real>(l) + 1.0) / (4.0 * M_PI)
             * _sigma_s_g_prime_g_l[_qp][scattering_index] * moment_l
             * _symmetry_factor;

      scattering_index++;
    }
  }

  return -1.0 * _test[_i][_qp] * res;
}
