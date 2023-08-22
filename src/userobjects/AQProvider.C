#include "AQProvider.h"

#include "GaussAngularQuadrature.h"

registerMooseObject("GnatApp", AQProvider);

InputParameters
AQProvider::validParams()
{
  auto params = ThreadedGeneralUserObject::validParams();
  params.addClassDescription("A user object which provides the angular quadrature set to discrete "
                             "ordinate transport schemes. This class should not be exposed to the "
                             "user, and should be enabled through a transport action.");
  params.addRequiredParam<MooseEnum>("dimensionality",
                                     MooseEnum("1D_cartesian 2D_cartesian 3D_cartesian"),
                                     "Dimensionality and the coordinate system of the "
                                     "problem.");
  params.addParam<MooseEnum>("aq_type",
                             MooseEnum("gauss_chebyshev", "gauss_chebyshev"),
                             "The angular quadrature set to use. Defaults to Gauss-Chebyshev.");
  params.addRangeCheckedParam<unsigned int>(
      "n_l",
      6,
      "n_l > 0",
      "Order of the polar product quadrature. Only required for Gauss-Chebyshev angular "
      "quadratures.");
  params.addRangeCheckedParam<unsigned int>(
      "n_c",
      6,
      "n_c > 0",
      "Order of the azimuthal product quadrature. Only required for Gauss-Chebyshev angular "
      "quadratures.");
  params.addParam<MooseEnum>("major_axis",
                             MooseEnum("x y z", "x"),
                             "Major axis of the angular quadrature. Allows the "
                             "polar angular quadrature to align with a cartesian "
                             "axis with minimal heterogeneity. Default is the "
                             "x-axis. This parameter is ignored for 1D and 2D "
                             "problems.");

  return params;
}

AQProvider::AQProvider(const InputParameters & parameters)
  : ThreadedGeneralUserObject(parameters),
    _aq_type(getParam<MooseEnum>("aq_type").getEnum<AQType>()),
    _aq(nullptr)
{
  switch (_aq_type)
  {
    case AQType::GaussChebyshev:
      _aq.reset(
          new GaussAngularQuadrature(getParam<unsigned int>("n_c"),
                                     getParam<unsigned int>("n_l"),
                                     getParam<MooseEnum>("major_axis").getEnum<MajorAxis>(),
                                     getParam<MooseEnum>("dimensionality").getEnum<ProblemType>()));
      break;
    default:
      break;
  }

  if (!_aq)
    mooseError("The angular quadrature set has not been initialized!");
}
