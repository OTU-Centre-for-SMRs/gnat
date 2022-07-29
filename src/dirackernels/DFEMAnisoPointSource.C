#include "DFEMAnisoPointSource.h"

#include "Function.h"

registerMooseObject("GnatApp", DFEMAnisoPointSource);

InputParameters
DFEMAnisoPointSource::validParams()
{
  auto params = SNBaseDiracKernel::validParams();
  params.addClassDescription("Computes the anisotropic point source term for the "
                             "current group of the discrete ordinates neutron "
                             "transport equation. The weak form is given by "
                             "$-(\\psi_{j}, S_{g}(\\hat{\\Omega}))$. "
                             "This kernel should not be exposed to the user, "
                             "instead being enabled through a transport action.");
  params.addRequiredRangeCheckedParam<unsigned int>("ordinate_index",
                                                    "ordinate_index >= 0",
                                                    "The discrete ordinate index "
                                                    "of the current angular "
                                                    "flux.");
  params.addRequiredParam<Real>("intensity",
                                "Intensity of the anisotropic point source.");
  params.addRequiredParam<FunctionName>("phase_function",
                                        "Anisotropic phase function for the "
                                        "point source. x, y, and z are treated "
                                        "as the components of the directional "
                                        "cosines along those axes.");
  params.addRequiredParam<Point>("point", "Location of the point source.");

  return params;
}

DFEMAnisoPointSource::DFEMAnisoPointSource(const InputParameters & parameters)
  : SNBaseDiracKernel(parameters)
  , _ordinate_index(getParam<unsigned int>("ordinate_index"))
  , _source_intensity(getParam<Real>("intensities"))
  , _angular_distribution(getFunction("phase_function"))
  , _source_location(getParam<Point>("points"))
{ }

void
DFEMAnisoPointSource::addPoints()
{
  addPoint(_source_location);
}

Real
DFEMAnisoPointSource::computeQpResidual()
{
  Point temp(_quadrature_set.direction(_ordinate_index)(0),
             _quadrature_set.direction(_ordinate_index)(1),
             _quadrature_set.direction(_ordinate_index)(2));

  // Hijacking the MOOSE function system so the user can parse in an analytical
  // phase function for this point source.
  return (-1.0 / M_PI) * _test[_i][_qp] * _angular_distribution.value(_t, temp)
         * _source_intensity * _symmetry_factor;
}

Real
DFEMAnisoPointSource::computeQpJacobian()
{
  return 0.0;
}
