#include "DFEMScattering.h"

#include "RealSphericalHarmonics.h"

registerMooseObject("GnatApp", DFEMScattering);

InputParameters
DFEMScattering::validParams()
{
  auto params = SNBaseKernel::validParams();
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

DFEMScattering::DFEMScattering(const InputParameters & parameters)
  : SNBaseKernel(parameters),
    _ordinate_index(getParam<unsigned int>("ordinate_index")),
    _group_index(getParam<unsigned int>("group_index")),
    _num_groups(getParam<unsigned int>("num_groups")),
    _max_anisotropy(getParam<unsigned int>("max_anisotropy")),
    _sigma_s_g_prime_g_l(getADMaterialProperty<std::vector<Real>>(
        getParam<std::string>("transport_system") + "scattering_matrix")),
    _anisotropy(getMaterialProperty<unsigned int>(getParam<std::string>("transport_system") +
                                                  "medium_anisotropy")),
    _num_dir_sh(0u)
{
  if (_group_index >= _num_groups)
    mooseError("The group index exceeds the number of energy groups.");

  if (_ordinate_index >= _quadrature_set.totalOrder())
    mooseError("The ordinates index exceeds the number of quadrature points.");

  const unsigned int num_coupled = coupledComponents("group_flux_ordinates");
  if (num_coupled != _quadrature_set.totalOrder() * _num_groups)
    mooseError("Mismatch between the angular flux ordinates and quadrature set.");

  // Fetch the flux ordinates and their derivatives.
  _group_flux_ordinate_ids.reserve(num_coupled);
  _group_flux_ordinates.reserve(num_coupled);
  for (unsigned int i = 0; i < num_coupled; ++i)
  {
    _group_flux_ordinate_ids.emplace_back(coupled("group_flux_ordinates", i));
    _group_flux_ordinates.emplace_back(&coupledValue("group_flux_ordinates", i));
  }

  // Pre-compute the spherical harmonics coefficients.
  // Handle different levels of dimensionality.
  Real mu = 0.0;
  Real omega = 0.0;
  switch (_quadrature_set.getProblemType())
  {
    // Legendre moments in 1D, looping over m is unecessary.
    case ProblemType::Cartesian1D:
      _num_dir_sh = _max_anisotropy + 1u;
      _y_l_m_n.reserve(_quadrature_set.totalOrder() * _num_dir_sh);

      for (unsigned int n = 0u; n < _quadrature_set.totalOrder(); ++n)
      {
        mu = _quadrature_set.getPolarRoot(n);
        omega = _quadrature_set.getAzimuthalAngularRoot(n);

        for (unsigned int l = 0u; l <= _max_anisotropy; ++l)
          _y_l_m_n.emplace_back(RealSphericalHarmonics::evaluate(l, 0, mu, omega));
      }
      break;

    // Need moments with m >= 0 for 2D.
    case ProblemType::Cartesian2D:
      _num_dir_sh = (_max_anisotropy + 1u) * (_max_anisotropy + 2u) / 2u;
      _y_l_m_n.reserve(_quadrature_set.totalOrder() * _num_dir_sh);

      for (unsigned int n = 0u; n < _quadrature_set.totalOrder(); ++n)
      {
        mu = _quadrature_set.getPolarRoot(n);
        omega = _quadrature_set.getAzimuthalAngularRoot(n);

        for (unsigned int l = 0u; l <= _max_anisotropy; ++l)
        {
          for (int m = 0; m <= static_cast<int>(l); ++m)
            _y_l_m_n.emplace_back(RealSphericalHarmonics::evaluate(l, m, mu, omega));
        }
      }
      break;

    // Need all moments in 3D.
    case ProblemType::Cartesian3D:
      _num_dir_sh = (_max_anisotropy + 1u) * (_max_anisotropy + 1u);
      _y_l_m_n.reserve(_quadrature_set.totalOrder() * _num_dir_sh);

      for (unsigned int n = 0u; n < _quadrature_set.totalOrder(); ++n)
      {
        mu = _quadrature_set.getPolarRoot(n);
        omega = _quadrature_set.getAzimuthalAngularRoot(n);

        for (unsigned int l = 0u; l <= _max_anisotropy; ++l)
        {
          for (int m = -1 * static_cast<int>(l); m <= static_cast<int>(l); ++m)
            _y_l_m_n.emplace_back(RealSphericalHarmonics::evaluate(l, m, mu, omega));
        }
      }
      break;

    default: // Defaults to doing nothing for now.
      break;
  }
}

Real
DFEMScattering::computeFluxMoment(unsigned int g_prime, unsigned int sh_offset)
{
  Real moment = 0.0;
  for (unsigned int n = 0; n < _quadrature_set.totalOrder(); ++n)
  {
    moment += _y_l_m_n[n * _num_dir_sh + sh_offset] *
              (*_group_flux_ordinates[g_prime * _quadrature_set.totalOrder() + n])[_qp] *
              _quadrature_set.weight(n);
  }

  return moment;
}

// Compute the full scattering term for both in-group and group-to-group
// scattering.
Real
DFEMScattering::computeQpResidual()
{
  // Quit early if no Legendre cross-section moments are provided.
  if (_sigma_s_g_prime_g_l[_qp].size() == 0u)
    return 0.0;

  // The maximum degree of anisotropy we can handle.
  const unsigned int max_anisotropy = std::min(_anisotropy[_qp], _max_anisotropy);

  unsigned int scattering_index = 0u;
  unsigned int sh_offset = 0u;

  Real res = 0.0;
  Real moment_l = 0.0;
  for (unsigned int g_prime = 0; g_prime < _num_groups; ++g_prime)
  {
    scattering_index = g_prime * _num_groups * _anisotropy[_qp] + _group_index * _anisotropy[_qp];

    for (unsigned int l = 0; l <= max_anisotropy; ++l)
    {
      // Handle different levels of dimensionality.
      switch (_quadrature_set.getProblemType())
      {
        // Legendre moments in 1D, looping over m is unecessary.
        case ProblemType::Cartesian1D:
          moment_l += computeFluxMoment(g_prime, sh_offset) *
                      _y_l_m_n[_ordinate_index * _num_dir_sh + sh_offset];
          sh_offset++;
          break;

        // Need moments with m >= 0 for 2D.
        case ProblemType::Cartesian2D:
          for (int m = 0; m <= static_cast<int>(l); ++m)
          {
            moment_l += computeFluxMoment(g_prime, sh_offset) *
                        _y_l_m_n[_ordinate_index * _num_dir_sh + sh_offset];
            sh_offset++;
          }
          break;

        // Need all moments in 3D.
        case ProblemType::Cartesian3D:
          for (int m = -1 * static_cast<int>(l); m <= static_cast<int>(l); ++m)
          {
            moment_l += computeFluxMoment(g_prime, sh_offset) *
                        _y_l_m_n[_ordinate_index * _num_dir_sh + sh_offset];
            sh_offset++;
          }
          break;

        default: // Defaults to doing nothing for now.
          break;
      }

      res += (2.0 * static_cast<Real>(l) + 1.0) / (4.0 * M_PI) *
             MetaPhysicL::raw_value(_sigma_s_g_prime_g_l[_qp][scattering_index + l]) * moment_l *
             _symmetry_factor;

      moment_l = 0.0;
    }

    sh_offset = 0u;
  }

  return -1.0 * _test[_i][_qp] * res;
}

Real
DFEMScattering::computeQpJacobian()
{
  // The maximum degree of anisotropy we can handle.
  const unsigned int max_anisotropy = std::min(_anisotropy[_qp], _max_anisotropy);
  // The current index into the scattering matrix.
  unsigned int scattering_index =
      _group_index * _num_groups * _anisotropy[_qp] + _group_index * _anisotropy[_qp];
  // The current index into the pre-computed SH functions.
  unsigned int sh_offset = 0u;

  Real jac = 0.0;
  Real d_moment_d_u = 0.0;

  // Handle different levels of dimensionality.
  switch (_quadrature_set.getProblemType())
  {
    case ProblemType::Cartesian1D:
      for (unsigned int l = 0u; l <= max_anisotropy; ++l)
      {
        d_moment_d_u += _y_l_m_n[_ordinate_index * _num_dir_sh + sh_offset] *
                        _y_l_m_n[_ordinate_index * _num_dir_sh + sh_offset] *
                        _quadrature_set.weight(_ordinate_index) * _phi[_j][_qp];

        jac += (2.0 * static_cast<Real>(l) + 1.0) / (4.0 * M_PI) *
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
          d_moment_d_u += _y_l_m_n[_ordinate_index * _num_dir_sh + sh_offset] *
                          _y_l_m_n[_ordinate_index * _num_dir_sh + sh_offset] *
                          _quadrature_set.weight(_ordinate_index) * _phi[_j][_qp];

          sh_offset++;
        }
        jac += (2.0 * static_cast<Real>(l) + 1.0) / (4.0 * M_PI) *
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
          d_moment_d_u += _y_l_m_n[_ordinate_index * _num_dir_sh + sh_offset] *
                          _y_l_m_n[_ordinate_index * _num_dir_sh + sh_offset] *
                          _quadrature_set.weight(_ordinate_index) * _phi[_j][_qp];

          sh_offset++;
        }
        jac += (2.0 * static_cast<Real>(l) + 1.0) / (4.0 * M_PI) *
               MetaPhysicL::raw_value(_sigma_s_g_prime_g_l[_qp][scattering_index + l]) *
               d_moment_d_u * _symmetry_factor;

        d_moment_d_u = 0.0;
      }
      break;

    default:
      break;
  }

  return -1.0 * _test[_i][_qp] * jac;
}

// TODO: Non-linear temperature and density need to be accounted for in the off-diagonal Jacobian.
Real
DFEMScattering::computeQpOffDiagJacobian(unsigned int jvar)
{
  unsigned int g_prime = 0u;
  unsigned int n_prime = 0u;
  bool not_found = true;

  for (unsigned int g = 0u; g < _num_groups; ++g)
  {
    for (unsigned int n = 0u; n < _quadrature_set.totalOrder(); ++n)
    {
      if (_group_flux_ordinate_ids[g * _num_groups + n] == jvar &&
          !(g == _group_index && n == _ordinate_index))
      {
        g_prime = g;
        n_prime = n;
        not_found = false;
      }
    }
  }

  if (not_found)
    return 0.0;

  // The maximum degree of anisotropy we can handle.
  const unsigned int max_anisotropy = std::min(_anisotropy[_qp], _max_anisotropy);
  // The current index into the scattering matrix.
  unsigned int scattering_index =
      g_prime * _num_groups * _anisotropy[_qp] + _group_index * _anisotropy[_qp];
  // The current index into the pre-computed SH functions.
  unsigned int sh_offset = 0u;

  Real jac = 0.0;
  Real d_moment_d_u = 0.0;

  // Handle different levels of dimensionality.
  switch (_quadrature_set.getProblemType())
  {
    case ProblemType::Cartesian1D:
      for (unsigned int l = 0u; l <= max_anisotropy; ++l)
      {
        d_moment_d_u += _y_l_m_n[_ordinate_index * _num_dir_sh + sh_offset] *
                        _y_l_m_n[n_prime * _num_dir_sh + sh_offset] *
                        _quadrature_set.weight(n_prime) * _phi[_j][_qp];

        jac += (2.0 * static_cast<Real>(l) + 1.0) / (4.0 * M_PI) *
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
          d_moment_d_u += _y_l_m_n[_ordinate_index * _num_dir_sh + sh_offset] *
                          _y_l_m_n[n_prime * _num_dir_sh + sh_offset] *
                          _quadrature_set.weight(n_prime) * _phi[_j][_qp];

          sh_offset++;
        }
        jac += (2.0 * static_cast<Real>(l) + 1.0) / (4.0 * M_PI) *
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
          d_moment_d_u += _y_l_m_n[_ordinate_index * _num_dir_sh + sh_offset] *
                          _y_l_m_n[n_prime * _num_dir_sh + sh_offset] *
                          _quadrature_set.weight(n_prime) * _phi[_j][_qp];

          sh_offset++;
        }
        jac += (2.0 * static_cast<Real>(l) + 1.0) / (4.0 * M_PI) *
               MetaPhysicL::raw_value(_sigma_s_g_prime_g_l[_qp][scattering_index + l]) *
               d_moment_d_u * _symmetry_factor;

        d_moment_d_u = 0.0;
      }
      break;

    default:
      break;
  }

  return -1.0 * _test[_i][_qp] * jac;
}
