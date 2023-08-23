#include "DiffusionVolumeSource.h"

registerMooseObject("GnatApp", DiffusionVolumeSource);

InputParameters
DiffusionVolumeSource::validParams()
{
  auto params = Kernel::validParams();
  params.addClassDescription(
      "Computes the source term for the "
      "transport equation simplified with the "
      "diffusion approximation, "
      "where the source moments are provided by the user. The weak form is given by "
      "$-(\\psi_{j}, S_{g})$. "
      "This kernel should not be exposed to the user, "
      "instead being enabled through a transport action.");
  params.addRequiredRangeCheckedParam<unsigned int>("group_index",
                                                    "group_index >= 0",
                                                    "The group index of the "
                                                    "current angular flux.");
  params.addRequiredRangeCheckedParam<unsigned int>("num_groups",
                                                    "num_groups >= 1",
                                                    "The number of spectral "
                                                    "energy groups.");
  params.addRequiredParam<std::vector<Real>>("group_source",
                                             "The external source moments for "
                                             "all energy groups.");

  return params;
}

DiffusionVolumeSource::DiffusionVolumeSource(const InputParameters & parameters)
  : Kernel(parameters),
    _group_index(getParam<unsigned int>("group_index")),
    _num_groups(getParam<unsigned int>("num_groups")),
    _source_moments(getParam<std::vector<Real>>("group_source"))
{
  if (_group_index >= _num_groups)
    mooseError("The group index exceeds the number of energy groups.");

  if (_source_moments.size() / _num_groups < 1u)
    mooseError("Not enough moments have been provided.");
}

Real
DiffusionVolumeSource::computeQpResidual()
{
  const unsigned int moment_index = _group_index * _source_moments.size() / _num_groups;

  return -1.0 * _test[_i][_qp] * _source_moments[moment_index];
}
