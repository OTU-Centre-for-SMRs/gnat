#include "FissionHeatSource.h"

registerMooseObject("GnatApp", FissionHeatSource);

InputParameters
FissionHeatSource::validParams()
{
  auto params = Kernel::validParams();
  params.addClassDescription(
      "A class which computes a fission heat production term for neutronics calculations coupled "
      "with heat conduction and/or thermo-hydraulics. The weak form is given by $-(\\phi, \\sum_{g "
      "= 1}^{G}\\kappa_{g}\\Sigma_{f,g}\\Phi_{g,0,0})$.");
  params.addRequiredRangeCheckedParam<unsigned int>(
      "num_groups",
      "num_groups > 0",
      "The number of spectral energy groups. This must match the number of provided fission "
      "heating cross-sections.");
  params.addRequiredCoupledVar("group_scalar_fluxes", "The scalar fluxes indexed by energy group.");
  params.addParam<std::string>(
      "transport_system",
      "",
      "Name of the transport system which will consume the provided material properties. If one is "
      "not provided the first transport system will be used.");
  params.addParam<Real>("scaling_factor", 1.0, "A scaling factor for the fission source.");

  return params;
}

FissionHeatSource::FissionHeatSource(const InputParameters & parameters)
  : Kernel(parameters),
    _num_groups(getParam<unsigned int>("num_groups")),
    _kappa_fission(getADMaterialProperty<std::vector<Real>>(
        getParam<std::string>("transport_system") + "heating_xs_g")),
    _scaling_factor(getParam<Real>("scaling_factor"))
{
  const unsigned int comp = coupledComponents("group_scalar_fluxes");
  if (comp != _num_groups)
    mooseError("The number of provided scalar fluxes does not match the number of provided energy "
               "groups.");

  _group_scalar_fluxes.reserve(comp);
  for (unsigned int g = 0u; g < comp; ++g)
    _group_scalar_fluxes.emplace_back(&coupledValue("group_scalar_fluxes", g));
}

Real
FissionHeatSource::computeQpResidual()
{
  Real res = 0.0;
  for (unsigned int g = 0u; g < _num_groups; ++g)
    res -= MetaPhysicL::raw_value(_kappa_fission[_qp][g]) * (*_group_scalar_fluxes[g])[_qp];

  return _scaling_factor * res;
}
