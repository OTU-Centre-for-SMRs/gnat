#include "SNSourceBC.h"

#include "RealSphericalHarmonics.h"

registerMooseObject("GnatApp", SNSourceBC);

InputParameters
SNSourceBC::validParams()
{
  auto params = SNBaseBC::validParams();
  params.addClassDescription("Computes the surface source boundary condition "
                             "with the weak form given by "
                             "$\\langle \\psi_{j},\\, \\hat{n}\\cdot"
                             "\\hat{\\Omega}\\Psi_{inc,\\, g}\\rangle_{\\Gamma_{i}}$, "
                             "$\\hat{n}\\cdot\\hat{\\Omega} \\leq 0$. "
                             "This kernel should not be exposed to the user, "
                             "instead being enabled through a transport action.");
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
  params.addRequiredParam<std::vector<Real>>("group_source",
                                             "The external source moments for "
                                             "all energy groups.");
  params.addParam<unsigned int>(
      "source_anisotropy", 0u, "The external source anisotropy of the medium.");

  return params;
}

SNSourceBC::SNSourceBC(const InputParameters & parameters)
  : SNBaseBC(parameters),
    _num_groups(getParam<unsigned int>("num_groups")),
    _group_index(getParam<unsigned int>("group_index")),
    _ordinate_index(getParam<unsigned int>("ordinate_index")),
    _source_moments(getParam<std::vector<Real>>("group_source")),
    _source_anisotropy(getParam<unsigned int>("source_anisotropy"))
{
  if (_group_index >= _num_groups)
    mooseError("The group index exceeds the number of energy groups.");

  if (_ordinate_index >= _aq.totalOrder())
    mooseError("The ordinates index exceeds the number of quadrature points.");

  switch (_mesh.dimension())
  {
    case 1u:
      _max_source_moments = (_source_anisotropy + 1u);
      _max_source_moments *= _num_groups;
      break;

    case 2u:
      _max_source_moments = (_source_anisotropy + 1u) * (_source_anisotropy + 2u) / 2u;
      _max_source_moments *= _num_groups;
      break;

    case 3u:
      _max_source_moments = (_source_anisotropy + 1u) * (_source_anisotropy + 1u);
      _max_source_moments *= _num_groups;
      break;

    default:
      mooseError("Unknown mesh dimensionality.");
      break;
  }

  // Warn the user if more parameters have been provided than required.
  if (_source_moments.size() > _max_source_moments)
  {
    mooseWarning("More source moments have been provided than possibly "
                 "supported with the given maximum source anisotropy and "
                 "number of groups. The vector will be truncated.");
  }

  // Error if the user did not provide enough parameters.
  if (_source_moments.size() < _max_source_moments)
    mooseError("Not enough source moments have been provided.");
}

Real
SNSourceBC::computeQpResidual()
{
  Real res = 0.0;
  Real src_l = 0.0;

  const Real & mu = _aq.getPolarRoot(_ordinate_index);
  const Real & omega = _aq.getAzimuthalAngularRoot(_ordinate_index);

  const unsigned int num_group_moments = _source_moments.size() / _num_groups;
  unsigned int moment_index = _group_index * num_group_moments;

  Real n_dot_omega = _aq.direction(_ordinate_index) * _normals[_qp];
  if (n_dot_omega >= 0.0)
    res += _u[_qp];
  else
  {
    for (unsigned int l = 0u; l <= _source_anisotropy; ++l)
    {
      // Handle different levels of dimensionality.
      switch (_aq.getProblemType())
      {
        // Legendre moments in 1D, looping over m is unecessary.
        case ProblemType::Cartesian1D:
          src_l +=
              _source_moments[moment_index] * RealSphericalHarmonics::evaluate(l, 0, mu, omega);
          moment_index++;
          break;

        // Need moments with m >= 0 for 2D.
        case ProblemType::Cartesian2D:
          for (int m = 0; m <= static_cast<int>(l); ++m)
          {
            src_l +=
                _source_moments[moment_index] * RealSphericalHarmonics::evaluate(l, m, mu, omega);
            moment_index++;
          }
          break;

        // Need all moments in 3D.
        case ProblemType::Cartesian3D:
          for (int m = -1 * static_cast<int>(l); m <= static_cast<int>(l); ++m)
          {
            src_l +=
                _source_moments[moment_index] * RealSphericalHarmonics::evaluate(l, m, mu, omega);
            moment_index++;
          }
          break;

        default: // Defaults to doing nothing for now.
          break;
      }

      res += src_l * (2.0 * static_cast<Real>(l) + 1.0) / (4.0 * libMesh::pi) * _symmetry_factor;
      src_l = 0.0;
    }
  }

  return res * n_dot_omega * _test[_i][_qp];
}

Real
SNSourceBC::computeQpJacobian()
{
  const auto n_dot_omega = _aq.direction(_ordinate_index) * _normals[_qp];
  return n_dot_omega >= 0.0 ? n_dot_omega * _phi[_j][_qp] * _test[_i][_qp] : 0.0;
}
