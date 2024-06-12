#include "Analytical1DSlabFlux.h"

registerMooseObject("GnatApp", Analytical1DSlabFlux);

InputParameters
Analytical1DSlabFlux ::validParams()
{
  auto params = Function::validParams();
  params.addClassDescription("A function which computes the analytical scalar flux for an "
                             "isotropic planar source embedded in a 1D infinite absorbing slab.");
  params.addParam<Real>("src_strength", 1.0, "The emission rate of the planar source.");
  params.addParam<Real>("src_location", 0.0, "The location of the planar source.");
  params.addParam<Real>("cross_section", 0.0, "The macroscopic total cross-section of the slab.");

  return params;
}

Analytical1DSlabFlux::Analytical1DSlabFlux(const InputParameters & parameters)
  : Function(parameters),
    _source_strength(getParam<Real>("src_strength")),
    _source_location(getParam<Real>("src_location")),
    _total_xs(getParam<Real>("cross_section"))
{
}

Real
Analytical1DSlabFlux::value(Real t, const Point & p) const
{
  constexpr Real euler_masch = 0.577215664901532860606512090082402431042159335939923598805;
  const Real tau = std::abs(p(0) - _source_location) * _total_xs;
  Real e_1 = -std::expint(-tau);

  return 0.5 * _source_strength * e_1;
}
