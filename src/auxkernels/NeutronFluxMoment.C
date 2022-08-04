#include "NeutronFluxMoment.h"

#include "RealSphericalHarmonics.h"

registerMooseObject("GnatApp", NeutronFluxMoment);

InputParameters
NeutronFluxMoment::validParams()
{
  auto params = AuxKernel::validParams();
  params.addClassDescription("Computes the flux moments "
                             "$\\Phi_{g,l,m}(\\vec{r}, t)$ of the scattering "
                             "source using the quadrature rule provided by "
                             "the material system.");
  params.addRequiredRangeCheckedParam<unsigned int>("n_l",
                                                    "n_l > 0",
                                                    "Order of the polar Gauss-"
                                                    "Legendre quadrature set.");
  params.addRequiredRangeCheckedParam<unsigned int>("n_c",
                                                    "n_c > 0",
                                                    "Order of the azimuthal "
                                                    "Gauss-Chebyshev "
                                                    "quadrature set.");
  params.addParam<MooseEnum>("major_axis", MooseEnum("x y z", "x"),
                             "Major axis of the angular quadrature. Allows the "
                             "polar angular quadrature to align with a cartesian "
                             "axis with minimal heterogeneity. Default is the "
                             "x-axis. This parameter is ignored for 1D and 2D "
                             "problems.");
  params.addRequiredParam<MooseEnum>("dimensionality",
                                     MooseEnum("1D_cartesian 2D_cartesian 3D_cartesian"),
                                     "Dimensionality and the coordinate system of the "
                                     "problem.");
  params.addRequiredCoupledVar("group_flux_ordinates", "The flux solutions for "
                               "all discrete directions. Must be listed in the "
                               "same order as the quadrature directions and "
                               "weights.");
  params.addRequiredParam<unsigned int>("degree", "Degree of this angular flux "
                                        "moment.");
  params.addRequiredParam<int>("order", "Order of this angular flux moment.");
  params.addParam<bool>("normalize_output", false,
                        "Divide the flux moment by the spherical harmonics "
                        "weights.");

  return params;
}

NeutronFluxMoment::NeutronFluxMoment(const InputParameters & parameters)
  : AuxKernel(parameters)
  , _quadrature_set(getParam<unsigned int>("n_c"),
                    getParam<unsigned int>("n_l"),
                    getParam<MooseEnum>("major_axis").getEnum<MajorAxis>(),
                    getParam<MooseEnum>("dimensionality").getEnum<ProblemType>())
  , _degree(getParam<unsigned int>("degree"))
  , _order(getParam<int>("order"))
  , _symmetry_factor(1.0)
{
  const unsigned int num_coupled = coupledComponents("group_flux_ordinates");

  if (num_coupled != _quadrature_set.totalOrder())
    mooseError("Mismatch between the angular flux ordinates and quadrature set.");

  _flux_ordinates.reserve(num_coupled);
  for (unsigned int i = 0; i < num_coupled; ++i)
    _flux_ordinates.emplace_back(&adCoupledValue("group_flux_ordinates", i));

  switch (_quadrature_set.getProblemType())
  {
    case ProblemType::Cartesian1D:
      _symmetry_factor = 2.0;
      break;

    case ProblemType::Cartesian2D:
      _symmetry_factor = 2.0 * M_PI;
      break;

    case ProblemType::Cartesian3D:
      _symmetry_factor = 4.0 * M_PI;
      break;

    default:
      _symmetry_factor = 1.0;
      break;
  }
}

void
NeutronFluxMoment::cartesianToSpherical(const RealVectorValue & ordinate,
                                        Real & mu, Real & omega)
{
  switch (_quadrature_set.getAxis())
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

Real
NeutronFluxMoment::computeValue()
{
  Real moment = 0.0;
  Real omega = 0.0;
  Real mu = 0.0;
  for (unsigned int i = 0; i < _quadrature_set.totalOrder(); ++i)
  {
    cartesianToSpherical(_quadrature_set.direction(i), mu, omega);
    moment += RealSphericalHarmonics::evaluate(_degree, _order, mu, omega)
              * std::max(MetaPhysicL::raw_value((* _flux_ordinates[i])[_qp]), 0.0)
              * _quadrature_set.weight(i);
  }

  return moment;
}
