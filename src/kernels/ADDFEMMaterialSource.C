#include "ADDFEMMaterialSource.h"

#include "RealSphericalHarmonics.h"

registerMooseObject("GnatApp", ADDFEMMaterialSource);

InputParameters
ADDFEMMaterialSource::validParams()
{
  auto params = ADSNBaseKernel::validParams();
  params.addClassDescription("Computes the source term for the "
                             "discrete ordinates neutron transport equation, "
                             "where the source moments are provided by the "
                             "material system. The weak form is given by "
                             "$-(\\psi_{j}, \\sum_{l = 0}^{L_{sr}} "
                             "\\frac{2l + 1}{4\\pi}\\sum_{m = -1}^{l} "
                             "S_{g,l,m}Y_{l,m}(\\hat{\\Omega}_{n}))$. "
                             "This kernel should not be exposed to the user, "
                             "instead being enabled through a transport action.");
  params.addRequiredRangeCheckedParam<unsigned int>("ordinate_index",
                                                    "ordinate_index >= 0",
                                                    "The discrete ordinate index "
                                                    "of the current angular "
                                                    "flux.");
  params.addRequiredRangeCheckedParam<unsigned int>("group_index",
                                                    "group_index >= 0",
                                                    "The group index of the "
                                                    "current angular flux.");
  params.addRequiredRangeCheckedParam<unsigned int>("num_groups",
                                                    "num_groups >= 1",
                                                    "The number of spectral "
                                                    "energy groups.");

  return params;
}

ADDFEMMaterialSource::ADDFEMMaterialSource(const InputParameters & parameters)
  : ADSNBaseKernel(parameters),
    _source_moments(getADMaterialProperty<std::vector<Real>>("source_moments")),
    _anisotropy(getMaterialProperty<unsigned int>("medium_source_anisotropy")),
    _ordinate_index(getParam<unsigned int>("ordinate_index")),
    _group_index(getParam<unsigned int>("group_index")),
    _num_groups(getParam<unsigned int>("num_groups"))
{
  if (_group_index >= _num_groups)
    mooseError("The group index exceeds the number of energy groups.");
}

ADReal
ADDFEMMaterialSource::computeQpResidual()
{
  if (_ordinate_index >= _quadrature_set.totalOrder())
    mooseError("The ordinates index exceeds the number of quadrature points.");

  // Quit early if there are no provided source moments.
  if (_source_moments[_qp].size() == 0u)
    return 0.0;

  const unsigned int num_group_moments = _source_moments[_qp].size() / _num_groups;

  ADReal src_l = 0.0;
  ADReal res = 0.0;
  Real mu = 0.0;
  Real omega = 0.0;
  unsigned int moment_index = _group_index * num_group_moments;
  for (unsigned int l = 0u; l <= _anisotropy[_qp]; ++l)
  {
    // Handle different levels of dimensionality.
    switch (_quadrature_set.getProblemType())
    {
      // Legendre moments in 1D, looping over m is unecessary.
      case ProblemType::Cartesian1D:
        cartesianToSpherical(_quadrature_set.direction(_ordinate_index), mu, omega);
        src_l +=
            _source_moments[_qp][moment_index] * RealSphericalHarmonics::evaluate(l, 0, mu, omega);
        moment_index++;
        break;

      // Need moments with m >= 0 for 2D.
      case ProblemType::Cartesian2D:
        for (int m = 0; m <= static_cast<int>(l); ++m)
        {
          cartesianToSpherical(_quadrature_set.direction(_ordinate_index), mu, omega);
          src_l += _source_moments[_qp][moment_index] *
                   RealSphericalHarmonics::evaluate(l, m, mu, omega);
          moment_index++;
        }
        break;

      // Need all moments in 3D.
      case ProblemType::Cartesian3D:
        for (int m = -1 * static_cast<int>(l); m <= static_cast<int>(l); ++m)
        {
          cartesianToSpherical(_quadrature_set.direction(_ordinate_index), mu, omega);
          src_l += _source_moments[_qp][moment_index] *
                   RealSphericalHarmonics::evaluate(l, m, mu, omega);
          moment_index++;
        }
        break;

      default: // Defaults to doing nothing for now.
        break;
    }

    res += src_l * (2.0 * static_cast<Real>(l) + 1.0) / (4.0 * M_PI) * _symmetry_factor;
    src_l = 0.0;
  }

  return -1.0 * _test[_i][_qp] * res;
}
