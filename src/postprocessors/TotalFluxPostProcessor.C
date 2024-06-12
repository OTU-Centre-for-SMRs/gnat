#include "TotalFluxPostProcessor.h"

#include "metaphysicl/raw_type.h"

registerMooseObject("GnatApp", TotalFluxPostProcessor);

InputParameters
TotalFluxPostProcessor::validParams()
{
  auto params = ElementIntegralPostprocessor::validParams();
  params.addClassDescription("A post-processor that computes the total (energy and volume "
                             "integrated) scalar flux of the radiation "
                             "transport equation.");
  params.addRequiredCoupledVar(
      "group_scalar_fluxes",
      "The scalar fluxes (zero'th moments of the angular fluxes) for all spectral energy groups.");
  params.addRequiredRangeCheckedParam<unsigned int>("num_groups",
                                                    "num_groups >= 1",
                                                    "The number of spectral "
                                                    "energy groups.");

  return params;
}

TotalFluxPostProcessor::TotalFluxPostProcessor(const InputParameters & parameters)
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
TotalFluxPostProcessor::computeQpIntegral()
{
  Real val = 0.0;
  for (unsigned int g = 0u; g < _num_groups; ++g)
    val += (*(_group_scalar_fluxes[g]))[_qp];

  return val;
}
