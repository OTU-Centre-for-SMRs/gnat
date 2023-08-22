#include "DiffusionIsoPointSource.h"

registerMooseObject("GnatApp", DiffusionIsoPointSource);

InputParameters
DiffusionIsoPointSource::validParams()
{
  auto params = DiracKernel::validParams();
  params.addClassDescription("Computes the isotropic point source term for "
                             "current group of the neutron "
                             "transport equation simplified with the "
                             "diffusion approximation. The weak form is given by "
                             "$-(\\psi_{j}, S_{g})$. "
                             "This kernel should not be exposed to the user, "
                             "instead being enabled through a transport action.");
  params.addRequiredRangeCheckedParam<unsigned int>("num_groups",
                                                    "num_groups >= 1",
                                                    "The number of spectral "
                                                    "energy groups.");
  params.addRequiredRangeCheckedParam<unsigned int>("group_index",
                                                    "group_index >= 0",
                                                    "The energy group index "
                                                    "of the current angular "
                                                    "flux.");
  params.addRequiredParam<Point>("point", "The location of the anisotropic point source.");
  params.addRequiredParam<std::vector<Real>>("group_source",
                                             "The external source moments for "
                                             "all energy groups.");
  params.addParam<unsigned int>(
      "source_anisotropy", 0u, "The external source anisotropy of the medium.");

  return params;
}

DiffusionIsoPointSource::DiffusionIsoPointSource(const InputParameters & parameters)
  : DiracKernel(parameters),
    _num_groups(getParam<unsigned int>("num_groups")),
    _group_index(getParam<unsigned int>("group_index")),
    _source_moments(getParam<std::vector<Real>>("group_source")),
    _source_location(getParam<Point>("point")),
    _anisotropy(getParam<unsigned int>("source_anisotropy"))
{
  if (_group_index >= _num_groups)
    mooseError("The group index exceeds the number of energy groups.");
}

void
DiffusionIsoPointSource::addPoints()
{
  addPoint(_source_location);
}

Real
DiffusionIsoPointSource::computeQpResidual()
{
  const unsigned int moment_index = _group_index * _source_moments.size() / _num_groups;

  return -1.0 * _test[_i][_qp] * _source_moments[moment_index];
}
