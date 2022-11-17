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
  params.addRequiredParam<std::vector<Real>>("intensities",
                                             "Isotropic source intensities for "
                                             "all isotropic point sources. "
                                             "Number of provided intensities "
                                             "must match the number of points "
                                             "in 'points'.");
  params.addRequiredParam<std::vector<Point>>("points",
                                              "All isotropic point source "
                                              "locations. Number of provided "
                                              "points must match the number of "
                                              "provided intensities in "
                                              "'intensities'.");

  return params;
}

DiffusionIsoPointSource::DiffusionIsoPointSource(const InputParameters & parameters)
  : DiracKernel(parameters),
    _source_intensities(getParam<std::vector<Real>>("intensities")),
    _source_locations(getParam<std::vector<Point>>("points"))
{
  if (_source_intensities.size() != _source_locations.size())
  {
    mooseError("The number of point source locations does not match the number "
               "of source intensities.");
  }
}

void
DiffusionIsoPointSource::addPoints()
{
  _point_intensity_mapping.clear();
  for (unsigned int i = 0; i < _source_intensities.size(); ++i)
  {
    _point_intensity_mapping[_source_locations[i]] = i;
    addPoint(_source_locations[i]);
  }
}

Real
DiffusionIsoPointSource::computeQpResidual()
{
  Real res = -1.0 * _test[_i][_qp] * _source_intensities[_point_intensity_mapping[_current_point]];

  return res;
}

Real
DiffusionIsoPointSource::computeQpJacobian()
{
  return 0.0;
}
