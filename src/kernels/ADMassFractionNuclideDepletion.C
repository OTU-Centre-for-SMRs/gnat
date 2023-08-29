#include "ADMassFractionNuclideDepletion.h"

registerMooseObject("GnatApp", ADMassFractionNuclideDepletion);

InputParameters
ADMassFractionNuclideDepletion::validParams()
{
  auto params = ADIsotopeBase::validParams();
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

ADMassFractionNuclideDepletion::ADMassFractionNuclideDepletion(const InputParameters & parameters)
  : ADIsotopeBase(parameters),
    _num_groups(getParam<unsigned int>("num_groups")),
    _sigma_a_g(getParam<std::vector<Real>>("group_absorption"))
{
  const unsigned int num_coupled = coupledComponents("group_scalar_fluxes");
  if (num_coupled != _sigma_a_g.size() && num_coupled != _num_groups)
  {
    mooseError("Mismatch between the number of provided cross-sections and the "
               "scalar fluxes.");
  }

  _group_scalar_fluxes.reserve(num_coupled);
  for (unsigned int i = 0u; i < num_coupled; ++i)
    _group_scalar_fluxes.emplace_back(&adCoupledValue("group_scalar_fluxes", i));
}

ADReal
ADMassFractionNuclideDepletion::computeQpResidual()
{
  auto qp_args = Moose::ElemQpArg();
  qp_args.elem = _current_elem;
  qp_args.qp = _qp;
  qp_args.qrule = _qrule;
  qp_args.point = _q_point[_qp];

  ADReal res = 0.0;

  // Loop over all scalar fluxes and microscopic cross-sections to compute the
  // sink.
  for (unsigned int g = 0u; g < _num_groups; ++g)
    res += (*_group_scalar_fluxes[g])[_qp] * _sigma_a_g[g];

  return computeQpTests() * res * _u[_qp] * _density(qp_args, 0u);
}
