#include "BaseNeutronicsMaterial.h"

registerMooseObject("GnatApp", BaseNeutronicsMaterial);

InputParameters
BaseNeutronicsMaterial::validParams()
{
  auto params = Material::validParams();
  params.addClassDescription("Provides basic functionality for the neutron "
                             "transport kernels, and should be applied over the "
                             "entire domain. This material handles the generation "
                             "of angular quadrature "
                             "sets. Note that these material properties should "
                             "not be exposed to the user, instead being enabled "
                             "through a transport action.");
  params.addRequiredRangeCheckedParam<unsigned int>("legendre_order",
                                                    "legendre_order > 0",
                                                    "Order of the polar Gauss-"
                                                    "Legendre quadrature set.");
  params.addRequiredRangeCheckedParam<unsigned int>("chebyshev_order",
                                                    "chebyshev_order > 0",
                                                    "Order of the azimuthal "
                                                    "Gauss-Chebyshev "
                                                    "quadrature set.");
  MooseEnum major_axis("x y z");
  params.addParam<MooseEnum>("major_axis", major_axis,
                             "Major axis of the angular quadrature. Allows the "
                             "polar angular quadrature to align with a cartesian "
                             "axis with minimal heterogeneity. Default is the "
                             "x-axis.");

  return params;
}

BaseNeutronicsMaterial::BaseNeutronicsMaterial(const InputParameters & parameters)
  : Material(parameters)
  , _quadrature_set(getParam<unsigned int>("chebyshev_order"),
                    getParam<unsigned int>("legendre_order"),
                    getParam<MooseEnum>("major_axis").getEnum<MajorAxis>())
  , _quadrature_directions(declareADProperty<std::vector<RealVectorValue>>("directions"))
  , _quadrature_weights(declareADProperty<std::vector<Real>>("direction_weights"))
  , _axis(declareProperty<MajorAxis>("quadrature_axis_alignment"))
{ }

void
BaseNeutronicsMaterial::computeQpProperties()
{
  auto & dir = _quadrature_set.getDirections();
  auto & w = _quadrature_set.getWeights();
  auto order = _quadrature_set.totalOrder();

  _quadrature_directions[_qp].resize(order, 0.0);
  _quadrature_weights[_qp].resize(order, 0.0);
  for (unsigned int i = 0; i < order; ++i)
  {
    _quadrature_directions[_qp][i] = dir[i];
    _quadrature_weights[_qp][i] = w[i];
  }

  _axis[_qp] = _quadrature_set.getAxis();
}
