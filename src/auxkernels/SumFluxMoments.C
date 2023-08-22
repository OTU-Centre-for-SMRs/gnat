#include "SumFluxMoments.h"

registerMooseObject("GnatApp", SumFluxMoments);

InputParameters
SumFluxMoments::validParams()
{
  auto params = AuxKernel::validParams();
  params.addClassDescription("A class which sums the uncollided and collided flux moments to yield "
                             "the total moments of the angular flux.");
  params.addRequiredCoupledVar("collided_flux_moment", "The collided moment of the angular flux.");
  params.addRequiredCoupledVar("uncollided_flux_moment",
                               "The uncollided moment of the angular flux.");

  return params;
}

SumFluxMoments::SumFluxMoments(const InputParameters & parameters)
  : AuxKernel(parameters),
    _collided_moment(adCoupledValue("collided_flux_moment")),
    _uncollided_moment(adCoupledValue("uncollided_flux_moment"))
{
}

Real
SumFluxMoments::computeValue()
{
  return MetaPhysicL::raw_value(_collided_moment[_qp]) +
         MetaPhysicL::raw_value(_uncollided_moment[_qp]);
}
