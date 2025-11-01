#include "SAAFVolumeSource.h"

#include "RealSphericalHarmonics.h"

registerMooseObject("GnatApp", SAAFVolumeSource);

InputParameters
SAAFVolumeSource::validParams()
{
  auto params = SAAFBaseKernel::validParams();
  params.addClassDescription("Computes the source term for the SAAF "
                             "discrete ordinates particle transport equation, "
                             "where the source moments are defined over a block. "
                             "The weak form is given by "
                             "$-(\\psi_{j} + \\tau_{g}\\vec{\\nabla}\\psi_{j}"
                             "\\cdot\\hat{\\Omega}, \\sum_{l = 0}^{L_{sr}} "
                             "\\frac{2l + 1}{4\\pi}\\sum_{m = -1}^{l} "
                             "S_{g,l,m}Y_{l,m}(\\hat{\\Omega}_{n}))$. "
                             "This kernel should not be exposed to the user, "
                             "instead being enabled through a transport action.");
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

SAAFVolumeSource::SAAFVolumeSource(const InputParameters & parameters)
  : SAAFBaseKernel(parameters),
    _num_groups(getParam<unsigned int>("num_groups")),
    _source_moments(getParam<std::vector<Real>>("group_source")),
    _anisotropy(getParam<unsigned int>("source_anisotropy"))
{
  if (_group_index >= _num_groups)
    mooseError("The group index exceeds the number of energy groups.");

  if (_ordinate_index >= _aq.totalOrder())
    mooseError("The ordinates index exceeds the number of quadrature points.");

  if (_anisotropy > 3)
    mooseError("Maximum degree of anisotropy supported is at most order 3 at present!");

  // Handle different levels of dimensionality.
  switch (_aq.getProblemType())
  {
    // Legendre moments in 1D, looping over m is unecessary.
    case ProblemType::Cartesian1D:
      if (_source_moments.size() > (_anisotropy + 1u) * _num_groups)
        mooseWarning("More source moments have been provided than possibly "
                     "supported with the given maximum source anisotropy and "
                     "number of groups. The vector will be truncated.");

      if (_source_moments.size() < (_anisotropy + 1u) * _num_groups)
        mooseError("Not enough source moments have been provided.");

      break;

    // Need moments with m >= 0 for 2D.
    case ProblemType::Cartesian2D:
      if (_source_moments.size() > ((_anisotropy + 1u) * (_anisotropy + 2u) / 2u) * _num_groups)
        mooseWarning("More source moments have been provided than possibly "
                     "supported with the given maximum source anisotropy and "
                     "number of groups. The vector will be truncated.");

      if (_source_moments.size() < ((_anisotropy + 1u) * (_anisotropy + 2u) / 2u) * _num_groups)
        mooseError("Not enough source moments have been provided.");

      break;

    // Need all moments in 3D.
    case ProblemType::Cartesian3D:
      if (_source_moments.size() > ((_anisotropy + 1u) * (_anisotropy + 1u)) * _num_groups)
        mooseWarning("More source moments have been provided than possibly "
                     "supported with the given maximum source anisotropy and "
                     "number of groups. The vector will be truncated.");

      if (_source_moments.size() < ((_anisotropy + 1u) * (_anisotropy + 1u)) * _num_groups)
        mooseError("Not enough source moments have been provided.");

      break;

    default: // Defaults to doing nothing for now.
      break;
  }
}

Real
SAAFVolumeSource::computeQpResidual()
{
  unsigned int moment_index = _group_index * _source_moments.size() / _num_groups;

  // The current quadrature direction.
  const auto & dir = _aq.direction(_ordinate_index);

  Real res = 0.0;
  for (unsigned int l = 0u; l <= _anisotropy; ++l)
  {
    Real src_l = 0.0;

    // Handle different levels of dimensionality.
    switch (_aq.getProblemType())
    {
      // Legendre moments in 1D, looping over m is unecessary.
      case ProblemType::Cartesian1D:
        src_l += _source_moments[moment_index] *
                 RealSphericalHarmonics::evaluate(l, 0, dir(0), dir(1), dir(2));
        moment_index++;
        break;

      // Need moments with m >= 0 for 2D.
      case ProblemType::Cartesian2D:
        for (int m = 0; m <= static_cast<int>(l); ++m)
        {
          src_l += _source_moments[moment_index] *
                   RealSphericalHarmonics::evaluate(l, m, dir(0), dir(1), dir(2));
          moment_index++;
        }
        break;

      // Need all moments in 3D.
      case ProblemType::Cartesian3D:
        for (int m = -1 * static_cast<int>(l); m <= static_cast<int>(l); ++m)
        {
          src_l += _source_moments[moment_index] *
                   RealSphericalHarmonics::evaluate(l, m, dir(0), dir(1), dir(2));
          moment_index++;
        }
        break;

      default: // Defaults to doing nothing for now.
        break;
    }

    res += src_l * (2.0 * static_cast<Real>(l) + 1.0) / (4.0 * libMesh::pi) * _symmetry_factor;
  }

  return -1.0 * computeQpTests() * res;
}
