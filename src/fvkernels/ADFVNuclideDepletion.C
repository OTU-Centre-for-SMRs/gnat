#include "ADFVNuclideDepletion.h"

registerMooseObject("GnatApp", ADFVNuclideDepletion);

InputParameters
ADFVNuclideDepletion::validParams()
{
  auto params = FVElementalKernel::validParams();
  params.addClassDescription("Computes the neutron absorption sink for the "
                             "isotope scalar transport equation: "
                             "$( \\phi_{j}, \\sum_{g = 1}^{G}"
                             "\\sigma_{a,g,i}N_{i}\\Phi_{g} )$.");
  params.addRequiredParam<unsigned int>("num_groups",
                                        "The number of spectral neutron energy "
                                        "groups.");
  params.addRequiredParam<std::vector<Real>>("group_absorption",
                                             "The microscopic neutron "
                                             "absorption cross-section for "
                                             "each energy group. These "
                                             "cross-sections must be listed in "
                                             "decending order by group.");
  params.addRequiredCoupledVar("group_scalar_fluxes",
                               "The scalar neutron fluxes for each energy "
                               "group.");

  return params;
}

ADFVNuclideDepletion::ADFVNuclideDepletion(const InputParameters & parameters)
  : FVElementalKernel(parameters),
    _num_groups(getParam<unsigned int>("num_groups")),
    _sigma_a_g(getParam<std::vector<Real>>("group_absorption"))
{
  const unsigned int num_coupled_fluxes = coupledComponents("group_scalar_fluxes");

  if (num_coupled_fluxes != _sigma_a_g.size() && num_coupled_fluxes != _num_groups)
  {
    mooseError("Mismatch between the number of provided cross-sections and the "
               "scalar fluxes.");
  }

  _group_scalar_fluxes.reserve(num_coupled_fluxes);
  for (unsigned int i = 0u; i < num_coupled_fluxes; ++i)
    _group_scalar_fluxes.emplace_back(&adCoupledValue("group_scalar_fluxes", i));
}

ADReal
ADFVNuclideDepletion::computeQpResidual()
{
  auto elem_args = makeElemArg(_current_elem);

  ADReal res = 0.0;
  // Loop over all scalar fluxes and microscopic cross-sections to compute the
  // sink.
  for (unsigned int g = 0u; g < _num_groups; ++g)
    res += ((*_group_scalar_fluxes[g])[_qp]) * _sigma_a_g[g];

  return res * _u_functor(elem_args, 0u);
}
