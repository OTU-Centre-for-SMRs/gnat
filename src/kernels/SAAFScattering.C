#include "SAAFScattering.h"

#include "RealSphericalHarmonics.h"

registerMooseObject("GnatApp", SAAFScattering);

InputParameters
SAAFScattering::validParams()
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
                             "moments are computed by this kernel. This kernel "
                             "should not be exposed to the user, instead being "
                             "enabled through a transport action.");
  params.addRequiredCoupledVar("group_flux_ordinates",
                               "The angular flux ordinates for all groups.");
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

SAAFScattering::SAAFScattering(const InputParameters & parameters)
  : SAAFBaseKernel(parameters),
    _num_groups(getParam<unsigned int>("num_groups")),
    _max_anisotropy(getParam<unsigned int>("max_anisotropy")),
    _sigma_s_g_prime_g_l(getADMaterialProperty<std::vector<Real>>(
        getParam<std::string>("transport_system") + "scattering_matrix")),
    _anisotropy(getMaterialProperty<unsigned int>(getParam<std::string>("transport_system") +
                                                  "medium_anisotropy"))
{
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

  if (_group_index >= _num_groups)
    mooseError("The group index exceeds the number of energy groups.");

  if (_ordinate_index >= _aq.totalOrder())
    mooseError("The ordinates index exceeds the number of quadrature points.");

  const unsigned int num_ord = coupledComponents("group_flux_ordinates");
  if (num_ord != _aq.totalOrder() * _num_groups)
    mooseError("Mismatch between the angular flux ordinates and quadrature set.");

  // Fetch the flux ordinate derivatives.
  for (unsigned int i = 0; i < num_ord; ++i)
  {
    unsigned int g = i / _aq.totalOrder();
    unsigned int n = i - g * _aq.totalOrder();
    _jvar_map.emplace(coupled("group_flux_ordinates", i), std::make_pair(g, n));
  }

  const unsigned int num_mom = coupledComponents("group_flux_moments");
  if (num_mom != _num_moments_per_group * _num_groups)
    mooseError("Mismatch between the number of angular flux moments and the provided anisotropy / "
               "number of groups.");

  // Fetch the flux moments.
  _group_flux_moments.reserve(num_mom);
  for (unsigned int i = 0; i < num_mom; ++i)
    _group_flux_moments.emplace_back(&coupledValue("group_flux_moments", i));

  if (_max_anisotropy > 3)
    mooseError("Maximum degree of anisotropy supported is at most order 3 at present!");
}

// Compute the full scattering term for both in-group and group-to-group
// scattering.
Real
SAAFScattering::computeQpResidual()
{
  // Quit early if no Legendre cross-section moments are provided.
  if (_sigma_s_g_prime_g_l[_qp].size() == 0u)
    return 0.0;

  // The maximum degree of anisotropy we can handle.
  const unsigned int max_anisotropy = std::min(_anisotropy[_qp], _max_anisotropy);

  // The current quadrature direction.
  const auto & dir = _aq.direction(_ordinate_index);

  Real res = 0.0;
  for (unsigned int g_prime = 0; g_prime < _num_groups; ++g_prime)
  {
    Real moment_l = 0.0;
    unsigned int sh_offset = 0u;
    unsigned int scattering_index =
        g_prime * _num_groups * (_anisotropy[_qp] + 1u) + _group_index * (_anisotropy[_qp] + 1u);

    for (unsigned int l = 0; l <= max_anisotropy; ++l)
    {
      // Handle different levels of dimensionality.
      switch (_aq.getProblemType())
      {
        // Legendre moments in 1D, looping over m is unecessary.
        case ProblemType::Cartesian1D:
          moment_l += (*_group_flux_moments[g_prime * _num_moments_per_group + sh_offset])[_qp] *
                      RealSphericalHarmonics::evaluate(l, 0, dir(0), dir(1), dir(2));
          sh_offset++;
          break;

        // Need moments with m >= 0 for 2D.
        case ProblemType::Cartesian2D:
          for (int m = 0; m <= static_cast<int>(l); ++m)
          {
            moment_l += (*_group_flux_moments[g_prime * _num_moments_per_group + sh_offset])[_qp] *
                        RealSphericalHarmonics::evaluate(l, m, dir(0), dir(1), dir(2));
            sh_offset++;
          }
          break;

        // Need all moments in 3D.
        case ProblemType::Cartesian3D:
          for (int m = -1 * static_cast<int>(l); m <= static_cast<int>(l); ++m)
          {
            moment_l += (*_group_flux_moments[g_prime * _num_moments_per_group + sh_offset])[_qp] *
                        RealSphericalHarmonics::evaluate(l, m, dir(0), dir(1), dir(2));
            sh_offset++;
          }
          break;

        default: // Defaults to doing nothing for now.
          break;
      }

      res += (2.0 * static_cast<Real>(l) + 1.0) / (4.0 * libMesh::pi) *
             MetaPhysicL::raw_value(_sigma_s_g_prime_g_l[_qp][scattering_index + l]) * moment_l *
             _symmetry_factor;
    }
  }

  return -1.0 * computeQpTests() * res;
}

Real
SAAFScattering::computeQpJacobian()
{
  // Quit early if no Legendre cross-section moments are provided.
  if (_sigma_s_g_prime_g_l[_qp].size() == 0u)
    return 0.0;

  // The maximum degree of anisotropy we can handle.
  const unsigned int max_anisotropy = std::min(_anisotropy[_qp], _max_anisotropy);
  // The current index into the scattering matrix.
  const unsigned int scattering_index =
      _group_index * _num_groups * (_anisotropy[_qp] + 1u) + _group_index * (_anisotropy[_qp] + 1u);

  // The current quadrature direction.
  const auto & dir = _aq.direction(_ordinate_index);
  // The current quadrature weight.
  const auto & aq_w = _aq.weight(_ordinate_index);

  Real jac = 0.0;
  Real d_moment_d_u = 0.0;

  // Handle different levels of dimensionality.
  switch (_aq.getProblemType())
  {
    case ProblemType::Cartesian1D:
      for (unsigned int l = 0u; l <= max_anisotropy; ++l)
      {
        d_moment_d_u += RealSphericalHarmonics::evaluate(l, 0, dir(0), dir(1), dir(2)) *
                        RealSphericalHarmonics::evaluate(l, 0, dir(0), dir(1), dir(2)) *
                        aq_w * _phi[_j][_qp];

        jac += (2.0 * static_cast<Real>(l) + 1.0) / (4.0 * libMesh::pi) *
               MetaPhysicL::raw_value(_sigma_s_g_prime_g_l[_qp][scattering_index + l]) *
               d_moment_d_u * _symmetry_factor;

        d_moment_d_u = 0.0;
      }
      break;

    case ProblemType::Cartesian2D:
      for (unsigned int l = 0u; l <= max_anisotropy; ++l)
      {
        for (int m = 0; m <= static_cast<int>(l); ++m)
        {
          d_moment_d_u += RealSphericalHarmonics::evaluate(l, m, dir(0), dir(1), dir(2)) *
                          RealSphericalHarmonics::evaluate(l, m, dir(0), dir(1), dir(2)) *
                          aq_w * _phi[_j][_qp];
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
          d_moment_d_u += RealSphericalHarmonics::evaluate(l, m, dir(0), dir(1), dir(2)) *
                          RealSphericalHarmonics::evaluate(l, m, dir(0), dir(1), dir(2)) *
                          aq_w * _phi[_j][_qp];
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

// TODO: Non-linear temperature and density need to be accounted for in the off-diagonal Jacobian.
Real
SAAFScattering::computeQpOffDiagJacobian(unsigned int jvar)
{
  // Quit early if no Legendre cross-section moments are provided.
  if (_sigma_s_g_prime_g_l[_qp].size() == 0u)
    return 0.0;

  auto [g_prime, n_prime] = _jvar_map[jvar];
  if (g_prime == _group_index && n_prime == _ordinate_index)
    return 0.0;

  // The maximum degree of anisotropy we can handle.
  const unsigned int max_anisotropy = std::min(_anisotropy[_qp], _max_anisotropy);
  // The current index into the scattering matrix.
  const unsigned int scattering_index =
      g_prime * _num_groups * (_anisotropy[_qp] + 1u) + _group_index * (_anisotropy[_qp] + 1u);

  // The outgoing quadrature direction.
  const auto & o_dir = _aq.direction(_ordinate_index);
  // The incoming quadrature direction.
  const auto & i_dir = _aq.direction(n_prime);
  // The incoming quadrature weight.
  const auto & i_aq_w = _aq.weight(n_prime);

  Real jac = 0.0;
  Real d_moment_d_u = 0.0;

  // Handle different levels of dimensionality.
  switch (_aq.getProblemType())
  {
    case ProblemType::Cartesian1D:
      for (unsigned int l = 0u; l <= max_anisotropy; ++l)
      {
        d_moment_d_u += RealSphericalHarmonics::evaluate(l, 0, o_dir(0), o_dir(1), o_dir(2)) *
                        RealSphericalHarmonics::evaluate(l, 0, i_dir(0), i_dir(1), i_dir(2)) *
                        i_aq_w * _phi[_j][_qp];

        jac += (2.0 * static_cast<Real>(l) + 1.0) / (4.0 * libMesh::pi) *
               MetaPhysicL::raw_value(_sigma_s_g_prime_g_l[_qp][scattering_index + l]) *
               d_moment_d_u * _symmetry_factor;

        d_moment_d_u = 0.0;
      }
      break;

    case ProblemType::Cartesian2D:
      for (unsigned int l = 0u; l <= max_anisotropy; ++l)
      {
        for (int m = 0; m <= static_cast<int>(l); ++m)
        {
          d_moment_d_u += RealSphericalHarmonics::evaluate(l, m, o_dir(0), o_dir(1), o_dir(2)) *
                          RealSphericalHarmonics::evaluate(l, m, i_dir(0), i_dir(1), i_dir(2)) *
                          i_aq_w * _phi[_j][_qp];
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
          d_moment_d_u += RealSphericalHarmonics::evaluate(l, m, o_dir(0), o_dir(1), o_dir(2)) *
                          RealSphericalHarmonics::evaluate(l, m, i_dir(0), i_dir(1), i_dir(2)) *
                          i_aq_w * _phi[_j][_qp];
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
