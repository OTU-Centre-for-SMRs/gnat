#include "NeutronFluxMoment.h"

registerMooseObject("GnatApp", NeutronFluxMoment);

InputParameters
NeutronFluxMoment::validParams()
{
  auto params = AuxKernel::validParams();
  params.addClassDescription("Computes the flux moments "
                             "$\\Phi_{g,l,m}(\\vec{r}, t)$ of the scattering "
                             "source using the quadrature rule provided by "
                             "the material system.");
  params.addRequiredCoupledVar("group_flux_ordinates", "The flux solutions for "
                               "all discrete directions. Must be listed in the "
                               "same order as the quadrature directions and "
                               "weights.");
  MooseEnum major_axis("x y z", "x");
  params.addParam<MooseEnum>("major_axis", major_axis,
                             "Major axis of the angular quadrature. Allows the "
                             "polar angular quadrature to align with a cartesian "
                             "axis with minimal heterogeneity. Default is the "
                             "x-axis. Must be equal to the major axis specified "
                             "in the BaseNeutronicsMaterial.");
  params.addRequiredParam<unsigned int>("degree", "Degree of this angular flux "
                                        "moment.");
  params.addRequiredParam<unsigned int>("order", "Order of this angular flux "
                                        "moment.");

  return params;
}

NeutronFluxMoment::NeutronFluxMoment(const InputParameters & parameters)
  : AuxKernel(parameters)
  , _quadrature_directions(getADMaterialProperty<std::vector<RealVectorValue>>("directions"))
  , _quadrature_weights(getADMaterialProperty<std::vector<Real>>("direction_weights"))
  , _axis(getParam<MooseEnum>("major_axis").getEnum<GaussAngularQuadrature::MajorAxis>())
  , _degree(getParam<unsigned int>("degree"))
  , _order(getParam<unsigned int>("order"))
{
  unsigned int num_coupled = coupledComponents("group_flux_ordinates");

  _flux_ordinates.reserve(num_coupled);
  for (unsigned int i = 0; i < num_coupled; ++i)
    _flux_ordinates.emplace_back(&coupledValue("group_flux_ordinates", i));
}

void
NeutronFluxMoment::cartesianToSpherical(RealVectorValue ordinate, Real & mu,
                                        Real & omega)
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

Real
NeutronFluxMoment::computeValue()
{
  if (_quadrature_directions[_qp].size() != _flux_ordinates.size()
      || _quadrature_weights[_qp].size() != _flux_ordinates.size())
  {
    mooseError("The number of flux ordiantes does not match the number of "
               "quadrature directions and/or weights.");
  }

  Real moment, omega, mu = 0.0;
  for (unsigned int i = 0; i < _quadrature_directions[_qp].size(); ++i)
  {
    cartesianToSpherical(MetaPhysicL::raw_value(_quadrature_directions[_qp][i]),
                         mu, omega);

    moment += RealSphericalHarmonics::evaluate(_degree, _order, mu, omega)
              * MetaPhysicL::raw_value((* _flux_ordinates[i])[_qp]
                                       * _quadrature_weights[_qp][i]);
  }

  return moment;
}
