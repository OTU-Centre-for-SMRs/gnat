#include "ADNeutronMaterialSource.h"

#include "RealSphericalHarmonics.h"

registerMooseObject("GnatApp", ADNeutronMaterialSource);

InputParameters
ADNeutronMaterialSource::validParams()
{
  auto params = ADKernel::validParams();
  params.addClassDescription("Computes the source term for the "
                             "discrete ordinates neutron transport equation, "
                             "where the source moments are provided by the "
                             "material system. The weak form is given by "
                             "$-(\\psi_{j}, \\sum_{l = 0}^{L_{sr}} "
                             "\\frac{2l + 1}{4\\pi}\\sum_{m = -1}^{l} "
                             "S_{g,l,m}Y_{l,m}(\\hat{\\Omega}_{n}))$. "
                             "This kernel should not be exposed to the user, "
                             "instead being enabled through a transport action.");
  MooseEnum major_axis("x y z", "x");
  params.addParam<MooseEnum>("major_axis", major_axis,
                             "Major axis of the angular quadrature. Allows the "
                             "polar angular quadrature to align with a cartesian "
                             "axis with minimal heterogeneity. Default is the "
                             "x-axis. Must be equal to the major axis specified "
                             "in the BaseNeutronicsMaterial.");
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

ADNeutronMaterialSource::ADNeutronMaterialSource(const InputParameters & parameters)
  : ADKernel(parameters)
  , _ordinate_index(getParam<unsigned int>("ordinate_index"))
  , _group_index(getParam<unsigned int>("group_index"))
  , _num_groups(getParam<unsigned int>("num_groups"))
  , _axis(getParam<MooseEnum>("major_axis").getEnum<GaussAngularQuadrature::MajorAxis>())
  , _source_moments(getADMaterialProperty<std::vector<Real>>("source_moments"))
  , _directions(getADMaterialProperty<std::vector<RealVectorValue>>("directions"))
{ }

void
ADNeutronMaterialSource::cartesianToSpherical(const RealVectorValue & ordinate,
                                              Real & mu, Real & omega)
{
  switch (_axis)
  {
    case GaussAngularQuadrature::MajorAxis::X:
      mu = ordinate(0);
      omega = std::acos(ordinate(1) / std::sqrt(1.0 - (mu * mu)));

      break;

    case GaussAngularQuadrature::MajorAxis::Y:
      mu = ordinate(1);
      omega = std::acos(ordinate(2) / std::sqrt(1.0 - (mu * mu)));

      break;

    case GaussAngularQuadrature::MajorAxis::Z:
      mu = ordinate(2);
      omega = std::acos(ordinate(0) / std::sqrt(1.0 - (mu * mu)));

      break;
  }
}

/*
 * We assume that the vector of source moments is stored in order of group first,
 * then moment indices. As an example for 2 energy groups (G = 2) and a 2nd
 * order real spherical harmonics source expansion with L = 2 is given below.
 *
 * The source moments are indexed as S_{g, l, m}:
 * _source_moments[_qp][0] = S_{1, 0, 0}
 * _source_moments[_qp][1] = S_{1, 1, -1}
 * _source_moments[_qp][2] = S_{1, 1, 0}
 * _source_moments[_qp][3] = S_{1, 1, 1}
 * _source_moments[_qp][4] = S_{1, 1, 1}
 * _source_moments[_qp][5] = S_{1, 2, -2}
 * _source_moments[_qp][6] = S_{1, 2, -1}
 * _source_moments[_qp][7] = S_{1, 2, 0}
 * _source_moments[_qp][8] = S_{1, 2, 1}
 * _source_moments[_qp][9] = S_{1, 2, 2}
 * _source_moments[_qp][10] = S_{2, 0, 0}
 * _source_moments[_qp][11] = S_{2, 1, -1}
 * _source_moments[_qp][12] = S_{2, 1, 0}
 * _source_moments[_qp][13] = S_{2, 1, 1}
 * _source_moments[_qp][14] = S_{2, 1, 1}
 * _source_moments[_qp][15] = S_{2, 2, -2}
 * _source_moments[_qp][16] = S_{2, 2, -1}
 * _source_moments[_qp][17] = S_{2, 2, 0}
 * _source_moments[_qp][18] = S_{2, 2, 1}
 * _source_moments[_qp][19] = S_{2, 2, 2}
 *
 * The material providing the source moments is expected to format them
 * according to this arrangement.
*/

ADReal
ADNeutronMaterialSource::computeQpResidual()
{
  if (_ordinate_index >= _directions[_qp].size())
    mooseError("The ordinates index exceeds the number of quadrature points.");
  if (_group_index >= _num_groups)
    mooseError("The group index exceeds the number of energy groups.");

  const unsigned int num_group_moments = _source_moments[_qp].size() / _num_groups;
  const unsigned int degree = std::sqrt(num_group_moments) - 1u;
  const unsigned int group_offset = _group_index * degree;

  ADReal res, src_l = 0.0;
  Real omega, mu = 0.0;

  unsigned int moment_index = group_offset;
  for (unsigned int l = 0u; l <= degree; ++l)
  {
    for (int m = -1 * static_cast<int>(l); m <= static_cast<int>(l); ++m)
    {
      cartesianToSpherical(MetaPhysicL::raw_value(_directions[_qp][_ordinate_index]),
                           mu, omega);
      src_l += _source_moments[_qp][moment_index]
               * RealSphericalHarmonics::evaluate(l, m, mu, omega);
      moment_index++;
    }
    res += src_l * (2.0 * static_cast<Real>(l) + 1.0) / (4.0 * M_PI);
    src_l = 0.0;
  }

  return -1.0 * _test[_i][_qp] * res;
}
