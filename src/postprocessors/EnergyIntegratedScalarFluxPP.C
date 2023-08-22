#include "EnergyIntegratedScalarFluxPP.h"

registerMooseObject("GnatApp", EnergyIntegratedScalarFluxPP);

InputParameters
EnergyIntegratedScalarFluxPP::validParams()
{
  auto params = ElementIntegralPostprocessor::validParams();
  params.addClassDescription("A post-processor that computes the energy integrated scalar flux.");
  params.addRequiredCoupledVar(
      "group_scalar_fluxes",
      "The scalar fluxes (zero'th moments of the angular fluxes) for all spectral energy groups.");
  params.addRequiredRangeCheckedParam<unsigned int>("num_groups",
                                                    "num_groups >= 1",
                                                    "The number of spectral "
                                                    "energy groups.");

  return params;
}

EnergyIntegratedScalarFluxPP::EnergyIntegratedScalarFluxPP(const InputParameters & parameters)
  : ElementIntegralPostprocessor(parameters), _num_groups(getParam<unsigned int>("num_groups"))
{
  const unsigned int num_coupled = coupledComponents("group_scalar_fluxes");
  if (num_coupled != _num_groups)
    mooseError("Mismatch between the number of scalar fluxes and the number of groups.");

  _group_scalar_fluxes.reserve(num_coupled);
  for (unsigned int i = 0u; i < num_coupled; ++i)
    _group_scalar_fluxes.emplace_back(&coupledValue("group_scalar_fluxes", i));
}

Real
EnergyIntegratedScalarFluxPP::computeQpIntegral()
{
  Real val = 0.0;
  for (unsigned int g = 0u; g < _num_groups; ++g)
    val += (*(_group_scalar_fluxes[g]))[_qp];

  return val;
}
