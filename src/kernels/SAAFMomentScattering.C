#include "SAAFMomentScattering.h"

#include "RealSphericalHarmonics.h"

registerMooseObject("GnatApp", SAAFMomentScattering);

InputParameters
SAAFMomentScattering::validParams()
{
  auto params = SAAFBaseKernel::validParams();
  params.addClassDescription("Computes the scattering term for the "
                             "current group of the SAAF discrete ordinates neutron "
                             "transport equation. The weak form is given by "
                             "$-(\\psi_{j} + \\tau_{g}\\vec{\\nabla}\\psi_{j}"
                             "\\cdot\\hat{\\Omega}, \\sum_{g' = 1}^{G}"
                             "\\Sigma_{s,\\, g'\\rightarrow g}"
                             "\\sum_{l = 0}^{L}\\frac{2l + 1}{4\\pi} "
                             "f_{g'\\rightarrow g,\\, l}"
                             "\\sum_{m = -l}^{l}Y_{l,m}(\\hat{\\Omega}_{n})"
                             "\\Phi_{g',l,m})$. The group flux "
                             "moments must be provided to this kernel. This kernel "
                             "should not be exposed to the user, instead being "
                             "enabled through a transport action.");
  params.addRequiredCoupledVar("group_flux_moments", "The angular flux moments for all groups.");
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

SAAFMomentScattering::SAAFMomentScattering(const InputParameters & parameters)
  : SAAFBaseKernel(parameters),
    _num_groups(getParam<unsigned int>("num_groups")),
    _max_anisotropy(getParam<unsigned int>("max_anisotropy")),
    _num_moments_per_group(0u),
    _sigma_s_g_prime_g_l(getADMaterialProperty<std::vector<Real>>(
        getParam<std::string>("transport_system") + "scattering_matrix")),
    _anisotropy(getMaterialProperty<unsigned int>(getParam<std::string>("transport_system") +
                                                  "medium_anisotropy"))
{
  if (_group_index >= _num_groups)
    mooseError("The group index exceeds the number of energy groups.");

  if (_ordinate_index >= _aq.totalOrder())
    mooseError("The ordinates index exceeds the number of quadrature points.");

  switch (_aq.getProblemType())
  {
    case ProblemType::Cartesian1D:
      _num_moments_per_group = _max_anisotropy + 1u;
      break;

    case ProblemType::Cartesian2D:
      _num_moments_per_group = (_max_anisotropy + 1u) * (_max_anisotropy + 2u) / 2u;
      break;

    case ProblemType::Cartesian3D:
      _num_moments_per_group = (_max_anisotropy + 1u) * (_max_anisotropy + 1u);
      break;

    default:
      _num_moments_per_group = 0u;
      break;
  }

  const unsigned int num_coupled = coupledComponents("group_flux_moments");
  if (num_coupled != _num_moments_per_group * _num_groups)
    mooseError("Mismatch between the number of angular flux moments and the provided anisotropy / "
               "number of groups.");

  // Fetch the flux ordinates and their derivatives.
  _group_flux_moments.reserve(num_coupled);
  for (unsigned int i = 0; i < num_coupled; ++i)
    _group_flux_moments.emplace_back(&coupledValue("group_flux_moments", i));

  // Pre-compute the spherical harmonics coefficients.
  // Handle different levels of dimensionality.
  const Real & mu = _aq.getPolarRoot(_ordinate_index);
  const Real & omega = _aq.getAzimuthalAngularRoot(_ordinate_index);
  _y_l_m.reserve(_num_moments_per_group);
  switch (_aq.getProblemType())
  {
    // Legendre moments in 1D, looping over m is unecessary.
    case ProblemType::Cartesian1D:
      for (unsigned int l = 0u; l <= _max_anisotropy; ++l)
        _y_l_m.emplace_back(RealSphericalHarmonics::evaluate(l, 0, mu, omega));

      break;

    // Need moments with m >= 0 for 2D.
    case ProblemType::Cartesian2D:
      for (unsigned int l = 0u; l <= _max_anisotropy; ++l)
      {
        for (int m = 0; m <= static_cast<int>(l); ++m)
          _y_l_m.emplace_back(RealSphericalHarmonics::evaluate(l, m, mu, omega));
      }

      break;

    // Need all moments in 3D.
    case ProblemType::Cartesian3D:
      for (unsigned int l = 0u; l <= _max_anisotropy; ++l)
      {
        for (int m = -1 * static_cast<int>(l); m <= static_cast<int>(l); ++m)
          _y_l_m.emplace_back(RealSphericalHarmonics::evaluate(l, m, mu, omega));
      }

      break;

    default: // Defaults to doing nothing for now.
      break;
  }
}

Real
SAAFMomentScattering::computeQpResidual()
{
  // Quit early if no Legendre cross-section moments are provided.
  if (_sigma_s_g_prime_g_l[_qp].size() == 0u)
    return 0.0;

  // The maximum degree of anisotropy we can handle.
  const unsigned int max_anisotropy = std::min(_anisotropy[_qp], _max_anisotropy);
  // The current index into the scattering matrix.
  unsigned int scattering_index = 0u;
  // The current index into the pre-computed SH functions.
  unsigned int sh_offset = 0u;

  Real res = 0.0;
  Real moment_l = 0.0;
  for (unsigned int g_prime = 0; g_prime < _num_groups; ++g_prime)
  {
    scattering_index =
        g_prime * _num_groups * (_anisotropy[_qp] + 1u) + _group_index * (_anisotropy[_qp] + 1u);

    for (unsigned int l = 0; l <= max_anisotropy; ++l)
    {
      // Handle different levels of dimensionality.
      switch (_aq.getProblemType())
      {
        // Legendre moments in 1D, looping over m is unecessary.
        case ProblemType::Cartesian1D:
          moment_l += (*_group_flux_moments[g_prime * _num_moments_per_group + sh_offset])[_qp] *
                      _y_l_m[sh_offset];
          sh_offset++;
          break;

        // Need moments with m >= 0 for 2D.
        case ProblemType::Cartesian2D:
          for (int m = 0; m <= static_cast<int>(l); ++m)
          {
            moment_l += (*_group_flux_moments[g_prime * _num_moments_per_group + sh_offset])[_qp] *
                        _y_l_m[sh_offset];
            sh_offset++;
          }
          break;

        // Need all moments in 3D.
        case ProblemType::Cartesian3D:
          for (int m = -1 * static_cast<int>(l); m <= static_cast<int>(l); ++m)
          {
            moment_l += (*_group_flux_moments[g_prime * _num_moments_per_group + sh_offset])[_qp] *
                        _y_l_m[sh_offset];
            sh_offset++;
          }
          break;

        default: // Defaults to doing nothing for now.
          break;
      }

      res += (2.0 * static_cast<Real>(l) + 1.0) / (4.0 * libMesh::pi) *
             MetaPhysicL::raw_value(_sigma_s_g_prime_g_l[_qp][scattering_index + l]) * moment_l *
             _symmetry_factor;

      moment_l = 0.0;
    }

    sh_offset = 0u;
  }

  return -1.0 * computeQpTests() * res;
}

// Assemble the within direction, within group Jacobian contribution.
Real
SAAFMomentScattering::computeQpJacobian()
{
  // Quit early if no Legendre cross-section moments are provided.
  if (_sigma_s_g_prime_g_l[_qp].size() == 0u)
    return 0.0;

  // The maximum degree of anisotropy we can handle.
  const unsigned int max_anisotropy = std::min(_anisotropy[_qp], _max_anisotropy);
  // The current index into the scattering matrix.
  unsigned int scattering_index =
      _group_index * _num_groups * (_anisotropy[_qp] + 1u) + _group_index * (_anisotropy[_qp] + 1u);
  // The current index into the pre-computed SH functions.
  unsigned int sh_offset = 0u;

  Real jac = 0.0;
  Real d_moment_d_u = 0.0;

  // Handle different levels of dimensionality.
  switch (_aq.getProblemType())
  {
    case ProblemType::Cartesian1D:
      for (unsigned int l = 0u; l <= max_anisotropy; ++l)
      {
        d_moment_d_u +=
            _y_l_m[sh_offset] * _y_l_m[sh_offset] * _aq.weight(_ordinate_index) * _phi[_j][_qp];

        jac += (2.0 * static_cast<Real>(l) + 1.0) / (4.0 * libMesh::pi) *
               MetaPhysicL::raw_value(_sigma_s_g_prime_g_l[_qp][scattering_index + l]) *
               d_moment_d_u * _symmetry_factor;

        d_moment_d_u = 0.0;
        sh_offset++;
      }
      break;

    case ProblemType::Cartesian2D:
      for (unsigned int l = 0u; l <= max_anisotropy; ++l)
      {
        for (int m = 0; m <= static_cast<int>(l); ++m)
        {
          d_moment_d_u +=
              _y_l_m[sh_offset] * _y_l_m[sh_offset] * _aq.weight(_ordinate_index) * _phi[_j][_qp];

          sh_offset++;
        }
        jac += (2.0 * static_cast<Real>(l) + 1.0) / (4.0 * libMesh::pi) *
               MetaPhysicL::raw_value(_sigma_s_g_prime_g_l[_qp][scattering_index + l]) *
               d_moment_d_u * _symmetry_factor;

        d_moment_d_u = 0.0;
      }
      break;

    case ProblemType::Cartesian3D:
      for (unsigned int l = 0u; l <= max_anisotropy; ++l)
      {
        for (int m = -1 * static_cast<int>(l); m <= static_cast<int>(l); ++m)
        {
          d_moment_d_u +=
              _y_l_m[sh_offset] * _y_l_m[sh_offset] * _aq.weight(_ordinate_index) * _phi[_j][_qp];

          sh_offset++;
        }
        jac += (2.0 * static_cast<Real>(l) + 1.0) / (4.0 * libMesh::pi) *
               MetaPhysicL::raw_value(_sigma_s_g_prime_g_l[_qp][scattering_index + l]) *
               d_moment_d_u * _symmetry_factor;

        d_moment_d_u = 0.0;
      }
      break;

    default:
      break;
  }

  return -1.0 * computeQpTests() * jac;
}
