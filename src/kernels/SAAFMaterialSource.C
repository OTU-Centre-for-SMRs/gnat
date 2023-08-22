#include "SAAFMaterialSource.h"

#include "RealSphericalHarmonics.h"

registerMooseObject("GnatApp", SAAFMaterialSource);

InputParameters
SAAFMaterialSource::validParams()
{
  auto params = SAAFBaseKernel::validParams();
  params.addClassDescription("Computes the source term for the SAAF "
                             "discrete ordinates neutron transport equation, "
                             "where the source moments are provided by the "
                             "material system. The weak form is given by "
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

  return params;
}

SAAFMaterialSource::SAAFMaterialSource(const InputParameters & parameters)
  : SAAFBaseKernel(parameters),
    _num_groups(getParam<unsigned int>("num_groups")),
    _source_moments(getADMaterialProperty<std::vector<Real>>(
        getParam<std::string>("transport_system") + "source_moments")),
    _anisotropy(getMaterialProperty<unsigned int>(getParam<std::string>("transport_system") +
                                                  "medium_source_anisotropy"))
{
  if (_group_index >= _num_groups)
    mooseError("The group index exceeds the number of energy groups.");

  if (_ordinate_index >= _aq.totalOrder())
    mooseError("The ordinates index exceeds the number of quadrature points.");
}

Real
SAAFMaterialSource::computeQpResidual()
{
  // Quit early if there are no provided source moments.
  if (_source_moments[_qp].size() == 0u)
    return 0.0;

  const unsigned int num_group_moments = _source_moments[_qp].size() / _num_groups;

  Real src_l = 0.0;
  Real res = 0.0;
  unsigned int moment_index = _group_index * num_group_moments;
  for (unsigned int l = 0u; l <= _anisotropy[_qp]; ++l)
  {
    // Handle different levels of dimensionality.
    switch (_aq.getProblemType())
    {
      // Legendre moments in 1D, looping over m is unecessary.
      case ProblemType::Cartesian1D:
        src_l += MetaPhysicL::raw_value(_source_moments[_qp][moment_index]) *
                 RealSphericalHarmonics::evaluate(l,
                                                  0,
                                                  _aq.getPolarRoot(_ordinate_index),
                                                  _aq.getAzimuthalAngularRoot(_ordinate_index));
        moment_index++;
        break;

      // Need moments with m >= 0 for 2D.
      case ProblemType::Cartesian2D:
        for (int m = 0; m <= static_cast<int>(l); ++m)
        {
          src_l += MetaPhysicL::raw_value(_source_moments[_qp][moment_index]) *
                   RealSphericalHarmonics::evaluate(l,
                                                    m,
                                                    _aq.getPolarRoot(_ordinate_index),
                                                    _aq.getAzimuthalAngularRoot(_ordinate_index));
          moment_index++;
        }
        break;

      // Need all moments in 3D.
      case ProblemType::Cartesian3D:
        for (int m = -1 * static_cast<int>(l); m <= static_cast<int>(l); ++m)
        {
          src_l += MetaPhysicL::raw_value(_source_moments[_qp][moment_index]) *
                   RealSphericalHarmonics::evaluate(l,
                                                    m,
                                                    _aq.getPolarRoot(_ordinate_index),
                                                    _aq.getAzimuthalAngularRoot(_ordinate_index));
          moment_index++;
        }
        break;

      default: // Defaults to doing nothing for now.
        break;
    }

    res += src_l * (2.0 * static_cast<Real>(l) + 1.0) / (4.0 * libMesh::pi) * _symmetry_factor;
    src_l = 0.0;
  }

  return -1.0 * computeQpTests() * res;
}
