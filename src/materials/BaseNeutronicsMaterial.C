#include "BaseNeutronicsMaterial.h"

registerMooseObject("GnatApp", BaseNeutronicsMaterial);

InputParameters
BaseNeutronicsMaterial::validParams()
{
  InputParameters params = Material::validParams();
  params.addClassDescription("Provides basic functionality for other neutron"
                             " transport kernels. This handles the generation"
                             " of Gauss-Chebyshev angular quadrature"
                             " set and spherical harmonics expansion"
                             " coefficients.");
  params.addRequiredRangeCheckedParam<unsigned int>("legendre_order",
                                                    "legendre_order > 0",
                                                    "Order of the polar Gauss-"
                                                    "Legendre quadrature set.");
  params.addRequiredRangeCheckedParam<unsigned int>("chebyshev_order",
                                                    "chebyshev_order > 0",
                                                    "Order of the azimuthal "
                                                    "Gauss-Chebyshev "
                                                    "quadrature set.");
  params.addRequiredRangeCheckedParam<unsigned int>("spherical_harmonics_order",
                                                    "spherical_harmonics_order >= 0",
                                                    "The spherical harmonics "
                                                    "expansion order for the "
                                                    "scattering and external "
                                                    "sources.");

  return params;
}

BaseNeutronicsMaterial::BaseNeutronicsMaterial(const InputParameters & parameters)
  : Material(parameters)
  , _harmonics_order(getParam<unsigned int>("spherical_harmonics_order"))
  , _quadrature_set(getParam<unsigned int>("chebyshev_order"),
                    getParam<unsigned int>("legendre_order"))
  , _quadrature_directions(declareADProperty<std::vector<RealVectorValue>>("directions"))
  , _quadrature_weights(declareADProperty<std::vector<Real>>("direction_weights"))
{
  if (_harmonics_order >= _quadrature_set.totalOrder() - 1)
    mooseWarning("Cannot fully integrate the requested SH order with the"
                 " provided quadrature order(s).");
}

void
BaseNeutronicsMaterial::computeQpProperties()
{
  auto & dir = _quadrature_set.getDirections();
  auto & w = _quadrature_set.getWeights();
  for (unsigned int i = 0; i < _quadrature_set.totalOrder(); ++i)
  {
    _quadrature_directions[_qp][i] = dir[i];
    _quadrature_weights[_qp][i] = w[i];
  }
}
