#include "ParticleFluxMomentSASF.h"

#include "RealSphericalHarmonics.h"

registerMooseObject("GnatApp", ParticleFluxMomentSASF);

InputParameters
ParticleFluxMomentSASF::validParams()
{
  auto params = AuxKernel::validParams();
  params.addClassDescription("An auxkernel which computes higher degree and order angular moments "
                             "of the uncollided flux in the Self-Adjoint Scalar Flux approach.");

  params.addRequiredParam<Point>("source_location", "The location of the point source.");

  params.addRequiredCoupledVar("uncollided_scalar_flux", "The uncollided scalar flux.");

  params.addRequiredParam<unsigned int>("degree",
                                        "Degree of this angular flux "
                                        "moment.");
  params.addRequiredParam<int>("order", "Order of this angular flux moment.");

  return params;
}

ParticleFluxMomentSASF::ParticleFluxMomentSASF(const InputParameters & parameters)
  : AuxKernel(parameters),
    _source_location(getParam<Point>("source_location")),
    _degree(getParam<unsigned int>("degree")),
    _order(getParam<int>("order")),
    _uncollided_scalar_flux(coupledValue("uncollided_scalar_flux"))
{
  if (_degree > 3)
    mooseError("Maximum degree of anisotropy supported is at most order 3 at present!");
}

// This is broken.
Real
ParticleFluxMomentSASF::computeValue()
{
  const auto direction =
      (_q_point[_qp] - _source_location) / ((_q_point[_qp] - _source_location).norm());

  return RealSphericalHarmonics::evaluate(_degree, _order, direction(0), direction(1), direction(2)) *
         _uncollided_scalar_flux[_qp];
}
