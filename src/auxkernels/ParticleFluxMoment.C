#include "ParticleFluxMoment.h"

#include "RealSphericalHarmonics.h"

registerMooseObject("GnatApp", ParticleFluxMoment);

InputParameters
ParticleFluxMoment::validParams()
{
  auto params = AuxKernel::validParams();
  params.addClassDescription("Computes the flux moments "
                             "$\\Phi_{g,l,m}(\\vec{r}, t)$ of the scattering "
                             "source using the quadrature rule provided by "
                             "the material system.");
  params.addRequiredParam<UserObjectName>(
      "aq", "The name of the angular quadrature provider user object.");
  params.addRequiredCoupledVar("group_flux_ordinates",
                               "The flux solutions for "
                               "all discrete directions. Must be listed in the "
                               "same order as the quadrature directions and "
                               "weights.");
  params.addRequiredParam<unsigned int>("degree",
                                        "Degree of this angular flux "
                                        "moment.");
  params.addRequiredParam<int>("order", "Order of this angular flux moment.");

  params.addParam<Real>("scale_factor", 1.0, "A scaling factor to apply to the flux moments.");

  params.addRequiredParam<unsigned int>("group_index", "The current spectral energy group.");
  params.addRequiredParam<unsigned int>("num_groups", "The number of spectral energy groups.");
  params.addCoupledVar("uncollided_flux_moments",
                       "The uncollided flux moments. Currently only supports uncollided scalar "
                       "fluxes.");

  return params;
}

ParticleFluxMoment::ParticleFluxMoment(const InputParameters & parameters)
  : AuxKernel(parameters),
    _aq(getUserObject<AQProvider>("aq")),
    _degree(getParam<unsigned int>("degree")),
    _order(getParam<int>("order")),
    _scale_factor(getParam<Real>("scale_factor")),
    _uncollided_scalar_flux(nullptr),
    _group_index(getParam<unsigned int>("group_index")),
    _num_groups(getParam<unsigned int>("num_groups"))
{
  const unsigned int num_coupled = coupledComponents("group_flux_ordinates");

  if (num_coupled != _aq.totalOrder())
    mooseError("Mismatch between the angular flux ordinates and quadrature set.");

  _flux_ordinates.reserve(num_coupled);
  for (unsigned int i = 0; i < num_coupled; ++i)
    _flux_ordinates.emplace_back(&adCoupledValue("group_flux_ordinates", i));

  if (isCoupled("uncollided_flux_moments"))
    _uncollided_scalar_flux = &coupledArrayValue("uncollided_flux_moments");
}

Real
ParticleFluxMoment::computeValue()
{
  Real moment = 0.0;

  // The collided component.
  for (unsigned int i = 0; i < _aq.totalOrder(); ++i)
  {
    moment += RealSphericalHarmonics::evaluate(
                  _degree, _order, _aq.getPolarRoot(i), _aq.getAzimuthalAngularRoot(i)) *
              MetaPhysicL::raw_value((*_flux_ordinates[i])[_qp]) * _aq.weight(i);
  }

  // The uncollided component.
  if (_uncollided_scalar_flux && _degree == 0u && _order == 0)
    moment += (*_uncollided_scalar_flux)[_qp](_group_index);

  return moment * _scale_factor;
}
