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
}

void
ParticleFluxMomentSASF::cartesianToSpherical(const RealVectorValue & direction,
                                             Real & mu,
                                             Real & omega)
{
  mu = direction(0);
  omega = 0.0;

  if (direction(1) > 0.0)
  {
    if (direction(2) > 0.0)
      omega += std::atan(std::abs(direction(2)) / std::abs(direction(1)));
    else if (direction(2) < 0.0)
      omega += std::atan(std::abs(direction(2)) / std::abs(direction(1))) + libMesh::pi;
  }
  else if (direction(1) < 0.0)
  {
    if (direction(2) > 0.0)
      omega += std::atan(std::abs(direction(2)) / std::abs(direction(1))) + 3.0 * libMesh::pi / 2.0;
    else if (direction(2) < 0.0)
      omega += std::atan(std::abs(direction(2)) / std::abs(direction(1))) + libMesh::pi / 2.0;
    else
      omega += libMesh::pi;
  }
  else
  {
    if (direction(2) > 0.0)
      omega += libMesh::pi / 2.0;
    else if (direction(2) < 0.0)
      omega += 3.0 * libMesh::pi / 2.0;
  }
}

// This is broken.
Real
ParticleFluxMomentSASF::computeValue()
{
  const auto direction =
      (_q_point[_qp] - _source_location) / ((_q_point[_qp] - _source_location).norm());

  Real mu = 0.0;
  Real omega = 0.0;
  cartesianToSpherical(direction, mu, omega);

  return RealSphericalHarmonics::evaluate(_degree, _order, mu, omega) *
         _uncollided_scalar_flux[_qp];
}
