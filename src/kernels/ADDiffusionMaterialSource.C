#include "ADDiffusionMaterialSource.h"

registerMooseObject("GnatApp", ADDiffusionMaterialSource);

InputParameters
ADDiffusionMaterialSource::validParams()
{
  auto params = ADKernel::validParams();
  params.addClassDescription("Computes the source term for the "
                             "transport equation simplified with the "
                             "diffusion approximation, "
                             "where the source moments are provided by the "
                             "material system. The weak form is given by "
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

  return params;
}

ADDiffusionMaterialSource::ADDiffusionMaterialSource(const InputParameters & parameters)
  : ADKernel(parameters),
    _group_index(getParam<unsigned int>("group_index")),
    _num_groups(getParam<unsigned int>("num_groups")),
    _source_moments(getADMaterialProperty<std::vector<Real>>("source_moments")),
    _anisotropy(getMaterialProperty<unsigned int>("medium_source_anisotropy"))
{
  if (_group_index >= _num_groups)
    mooseError("The group index exceeds the number of energy groups.");
}

ADReal
ADDiffusionMaterialSource::computeQpResidual()
{
  // Quit early if there are no provided source moments.
  if (_source_moments[_qp].size() == 0u)
    return 0.0;

  const unsigned int num_group_moments = _source_moments[_qp].size() / _num_groups;
  return -1.0 * _test[_i][_qp] * _source_moments[_qp][_group_index * num_group_moments];
}
