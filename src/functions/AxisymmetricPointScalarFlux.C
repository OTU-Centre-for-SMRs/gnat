#include "AxisymmetricPointScalarFlux.h"

registerMooseObject("GnatApp", AxisymmetricPointScalarFlux);

InputParameters
AxisymmetricPointScalarFlux::validParams()
{
  auto params = Function::validParams();
  params.addClassDescription("A function which computes the analytical scalar flux for a point "
                             "source surrounded by axisymmetric shielding media.");
  params.addRequiredParam<unsigned int>(
      "num_dims",
      "The dimensionality of the solution. 2D yields a flatland, 3D yields spherical cartesian.");
  params.addParam<Real>("src_strength", 1.0, "The emission rate of the point source.");
  params.addParam<Point>("src_location", Point(0.0, 0.0, 0.0), "The location of the point source.");
  params.addParam<std::vector<Real>>(
      "axisymmetric_radii", std::vector<Real>(), "The ending radii for each axisymmetric region.");
  params.addParam<std::vector<Real>>(
      "axisymmetric_cross_sections",
      std::vector<Real>(),
      "The macroscopic total cross-sections for each radial region.");

  return params;
}

AxisymmetricPointScalarFlux::AxisymmetricPointScalarFlux(const InputParameters & parameters)
  : Function(parameters),
    _dims(getParam<unsigned int>("num_dims")),
    _source_strength(getParam<Real>("src_strength")),
    _source_location(getParam<Point>("src_location")),
    _radii(getParam<std::vector<Real>>("axisymmetric_radii")),
    _total_xs(getParam<std::vector<Real>>("axisymmetric_cross_sections"))
{
  if (_total_xs.size() != _radii.size())
    mooseError(
        "Mismatch between the number of provided radii and the number of provided cross-sections.");
}

Real
AxisymmetricPointScalarFlux::value(Real t, const Point & p) const
{
  const Real r = (p - _source_location).norm();

  // Streaming component.
  Real val =
      _source_strength / (_dims == 3u ? (4.0 * libMesh::pi * r * r) : (2.0 * libMesh::pi * r));

  // Attenuation over each region.
  Real prev_radius = 0.0;
  for (unsigned int i = 0u; i < _radii.size(); ++i)
  {
    if (r > _radii[i])
    {
      val *= std::exp(-1.0 * _total_xs[i] * (_radii[i] - prev_radius));
      prev_radius = _radii[i];
    }
    else if (r <= _radii[i])
    {
      val *= std::exp(-1.0 * _total_xs[i] * (r - prev_radius));
      break;
    }
  }

  return val;
}
