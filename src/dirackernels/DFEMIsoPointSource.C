#include "DFEMIsoPointSource.h"

registerMooseObject("GnatApp", DFEMIsoPointSource);

InputParameters
DFEMIsoPointSource::validParams()
{
  auto params = SNBaseDiracKernel::validParams();
  params.addClassDescription("Computes the isotropic point source term for "
                             "current group of the discrete ordinates neutron "
                             "transport equation. The weak form is given by "
                             "$-(\\psi_{j}, \\frac{S_{g,0,0}}{4\\pi})$. "
                             "This kernel should not be exposed to the user, "
                             "instead being enabled through a transport action.");
  params.addRequiredParam<std::vector<Real>>("intensities",
                                             "Isotropic source intensities for "
                                             "all isotropic point sources. "
                                             "Number of provided intensities "
                                             "must match the number of points "
                                             "in 'points'.");
  params.addRequiredParam<std::vector<Point>>("points", "All isotropic point source "
                                              "locations. Number of provided "
                                              "points must match the number of "
                                              "provided intensities in "
                                              "'intensities'.");

  return params;
}

DFEMIsoPointSource::DFEMIsoPointSource(const InputParameters & parameters)
  : SNBaseDiracKernel(parameters)
  , _source_intensities(getParam<std::vector<Real>>("intensities"))
  , _source_locations(getParam<std::vector<Point>>("points"))
{
  if (_source_intensities.size() != _source_locations.size())
  {
    mooseError("The number of point source locations does not match the number "
               "of source intensities.");
  }
}

void
DFEMIsoPointSource::addPoints()
{
  _point_intensity_mapping.clear();
  for (unsigned int i = 0; i < _source_intensities.size(); ++i)
  {
    _point_intensity_mapping[_source_locations[i]] = i;
    addPoint(_source_locations[i]);
  }
}

Real
DFEMIsoPointSource::computeQpResidual()
{
  Real res = (-1.0 / (4.0 * M_PI)) * _test[_i][_qp]
             * _source_intensities[_point_intensity_mapping[_current_point]]
             * _symmetry_factor;

  return res;
}

Real
DFEMIsoPointSource::computeQpJacobian()
{
  return 0.0;
}
