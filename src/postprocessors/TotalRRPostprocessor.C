#include "TotalRRPostprocessor.h"

#include "metaphysicl/raw_type.h"

registerMooseObject("GnatApp", TotalRRPostprocessor);

InputParameters
TotalRRPostprocessor::validParams()
{
  auto params = ElementIntegralPostprocessor::validParams();
  params.addClassDescription(
      "A post-processor that computes the total reaction rate of the radiation "
      "transport equation.");
  params.addParam<std::string>(
      "transport_system",
      "",
      "Name of the transport system which will consume the provided material properties. If one is "
      "not provided the first transport system will be used.");
  params.addRequiredCoupledVar(
      "group_scalar_fluxes",
      "The scalar fluxes (zero'th moments of the angular fluxes) for all spectral energy groups.");
  params.addRequiredRangeCheckedParam<unsigned int>("num_groups",
                                                    "num_groups >= 1",
                                                    "The number of spectral "
                                                    "energy groups.");

  return params;
}

TotalRRPostprocessor::TotalRRPostprocessor(const InputParameters & parameters)
  : ElementIntegralPostprocessor(parameters),
    _num_groups(getParam<unsigned int>("num_groups")),
    _sigma_t_g(getADMaterialProperty<std::vector<Real>>(getParam<std::string>("transport_system") +
                                                        "total_xs_g"))
{
  const unsigned int num_coupled = coupledComponents("group_scalar_fluxes");
  if (num_coupled != _num_groups)
    mooseError("Mismatch between the number of scalar fluxes and the number of groups.");

  _group_scalar_fluxes.reserve(num_coupled);
  for (unsigned int i = 0u; i < num_coupled; ++i)
    _group_scalar_fluxes.emplace_back(&coupledValue("group_scalar_fluxes", i));
}

Real
TotalRRPostprocessor::computeQpIntegral()
{
  Real val = 0.0;
  for (unsigned int g = 0u; g < _num_groups; ++g)
    val += MetaPhysicL::raw_value(_sigma_t_g[_qp][g]) * (*(_group_scalar_fluxes[g]))[_qp];

  return val;
}
