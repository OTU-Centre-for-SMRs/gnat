#include "DiffusionMaterialSource.h"

registerMooseObject("GnatApp", DiffusionMaterialSource);

InputParameters
DiffusionMaterialSource::validParams()
{
  auto params = Kernel::validParams();
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
  params.addParam<std::string>(
      "transport_system",
      "",
      "Name of the transport system which will consume the provided material properties. If one is "
      "not provided the first transport system will be used.");

  return params;
}

DiffusionMaterialSource::DiffusionMaterialSource(const InputParameters & parameters)
  : Kernel(parameters),
    _group_index(getParam<unsigned int>("group_index")),
    _num_groups(getParam<unsigned int>("num_groups")),
    _source_moments(getADMaterialProperty<std::vector<Real>>(
        getParam<std::string>("transport_system") + "source_moments"))
{
  if (_group_index >= _num_groups)
    mooseError("The group index exceeds the number of energy groups.");
}

Real
DiffusionMaterialSource::computeQpResidual()
{
  // Quit early if there are no provided source moments.
  if (_source_moments[_qp].size() == 0u)
    return 0.0;

  const unsigned int num_group_moments = _source_moments[_qp].size() / _num_groups;
  return -1.0 * _test[_i][_qp] *
         MetaPhysicL::raw_value(_source_moments[_qp][_group_index * num_group_moments]);
}
