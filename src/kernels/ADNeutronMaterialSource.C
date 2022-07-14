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
  , _source_moments(getADMaterialProperty<std::vector<Real>>("source_moments"))
  , _directions(getADMaterialProperty<std::vector<RealVectorValue>>("directions"))
  , _axis(getMaterialProperty<MajorAxis>("quadrature_axis_alignment"))
  , _anisotropy(getMaterialProperty<unsigned int>("medium_source_anisotropy"))
  , _ordinate_index(getParam<unsigned int>("ordinate_index"))
  , _group_index(getParam<unsigned int>("group_index"))
  , _num_groups(getParam<unsigned int>("num_groups"))
{
  if (_group_index >= _num_groups)
    mooseError("The group index exceeds the number of energy groups.");
}

void
ADNeutronMaterialSource::cartesianToSpherical(const RealVectorValue & ordinate,
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
ADNeutronMaterialSource::computeQpResidual()
{
  if (_ordinate_index >= _directions[_qp].size())
    mooseError("The ordinates index exceeds the number of quadrature points.");

  const unsigned int num_group_moments = _source_moments[_qp].size() / _num_groups;

  ADReal res, src_l = 0.0;
  Real omega, mu = 0.0;

  unsigned int moment_index = _group_index * num_group_moments;
  for (unsigned int l = 0u; l <= _anisotropy[_qp]; ++l)
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
