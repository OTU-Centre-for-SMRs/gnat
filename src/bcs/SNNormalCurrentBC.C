#include "SNNormalCurrentBC.h"

#include "RealSphericalHarmonics.h"

registerMooseObject("GnatApp", SNNormalCurrentBC);

InputParameters
SNNormalCurrentBC::validParams()
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
  params.addRequiredParam<Real>("group_current", "The currents for the group.");
  params.addParam<unsigned int>(
      "current_anisotropy", 0u, "The anisotropy of the external current.");

  return params;
}

SNNormalCurrentBC::SNNormalCurrentBC(const InputParameters & parameters)
  : SNBaseBC(parameters),
    _num_groups(getParam<unsigned int>("num_groups")),
    _group_index(getParam<unsigned int>("group_index")),
    _ordinate_index(getParam<unsigned int>("ordinate_index")),
    _current(getParam<Real>("group_current")),
    _current_anisotropy(getParam<unsigned int>("current_anisotropy"))
{
  if (_group_index >= _num_groups)
    mooseError("The group index exceeds the number of energy groups.");

  if (_ordinate_index >= _aq.totalOrder())
    mooseError("The ordinates index exceeds the number of quadrature points.");
}

// Use the sifting property of the Dirac delta function to get rid of the integral and evaluate it
// at the boundary normal vector.
Real
SNNormalCurrentBC::computeQpResidual()
{
  Real res = 0.0;
  Real src_l = 0.0;

  const Real & mu = _aq.getPolarRoot(_ordinate_index);
  const Real & omega = _aq.getAzimuthalAngularRoot(_ordinate_index);

  Real n_mu = 0.0;
  Real n_omega = 0.0;
  cartesianToSpherical(-1.0 * MetaPhysicL::raw_value(_normals[_qp]), n_mu, n_omega);

  Real n_dot_omega = _aq.direction(_ordinate_index) * _normals[_qp];
  if (n_dot_omega >= 0.0)
    res += _u[_qp];
  else
  {
    for (unsigned int l = 0u; l <= _current_anisotropy; ++l)
    {
      // Handle different levels of dimensionality.
      switch (_aq.getProblemType())
      {
        // Legendre moments in 1D, looping over m is unecessary.
        case ProblemType::Cartesian1D:
          src_l += RealSphericalHarmonics::evaluate(l, 0, mu, omega) * _current *
                   RealSphericalHarmonics::evaluate(l, 0, n_mu, n_omega);
          break;

        // Need moments with m >= 0 for 2D.
        case ProblemType::Cartesian2D:
          for (int m = 0; m <= static_cast<int>(l); ++m)
          {
            src_l += RealSphericalHarmonics::evaluate(l, m, mu, omega) * _current *
                     RealSphericalHarmonics::evaluate(l, m, n_mu, n_omega);
          }
          break;

        // Need all moments in 3D.
        case ProblemType::Cartesian3D:
          for (int m = -1 * static_cast<int>(l); m <= static_cast<int>(l); ++m)
          {
            src_l += RealSphericalHarmonics::evaluate(l, m, mu, omega) * _current *
                     RealSphericalHarmonics::evaluate(l, m, n_mu, n_omega);
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
SNNormalCurrentBC::computeQpJacobian()
{
  const auto n_dot_omega = _aq.direction(_ordinate_index) * _normals[_qp];
  return n_dot_omega >= 0.0 ? n_dot_omega * _phi[_j][_qp] * _test[_i][_qp] : 0.0;
}
