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
  , _directions(getADMaterialProperty<std::vector<RealVectorValue>>("directions"))
  , _axis(getMaterialProperty<MajorAxis>("quadrature_axis_alignment"))
  , _sigma_s_g_prime_g_l(getADMaterialProperty<std::vector<Real>>("scattering_matrix"))
  , _anisotropy(getMaterialProperty<unsigned int>("medium_anisotropy"))
  , _ordinate_index(getParam<unsigned int>("ordinate_index"))
  , _group_index(getParam<unsigned int>("group_index"))
  , _num_groups(getParam<unsigned int>("num_groups"))
{
  if (_group_index >= _num_groups)
    mooseError("The group index exceeds the number of energy groups.");

  const unsigned int num_coupled = coupledComponents("group_flux_moments");

  const unsigned int num_group_moments = num_coupled / _num_groups;
  _provided_moment_degree = std::sqrt(num_group_moments) - 1u;

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

  ADReal res, moment_l = 0.0;
  Real omega, mu = 0.0;

  unsigned int scattering_index = 0u;
  unsigned int moment_index = 0u;
  for (unsigned int g_prime = 0; g_prime < _num_groups; ++g_prime)
  {
    moment_index = g_prime * (_provided_moment_degree + 1u) * (_provided_moment_degree + 1u);
    scattering_index = g_prime * _num_groups * _anisotropy[_qp] + _group_index * _anisotropy[_qp];

    if (g_prime == _group_index)
      continue;

    for (unsigned int l = 0; l <= _anisotropy[_qp]; ++l)
    {
      for (int m = -1 * static_cast<int>(l); m <= static_cast<int>(l); ++m)
      {
        cartesianToSpherical(MetaPhysicL::raw_value(_directions[_qp][_ordinate_index]),
                             mu, omega);
        moment_l += (*_group_flux_moments[moment_index])[_qp]
                    * RealSphericalHarmonics::evaluate(l, m, mu, omega);
        moment_index++;
      }

      res += (2.0 * static_cast<Real>(l) + 1.0) / (4.0 * M_PI)
             * _sigma_s_g_prime_g_l[_qp][scattering_index] * moment_l;

      scattering_index++;
    }
  }

  return -1.0 * _test[_i][_qp] * res;
}
