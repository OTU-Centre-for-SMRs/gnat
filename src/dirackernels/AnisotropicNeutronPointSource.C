#include "AnisotropicNeutronPointSource.h"

#include "Function.h"

registerMooseObject("GnatApp", AnisotropicNeutronPointSource);

InputParameters
AnisotropicNeutronPointSource::validParams()
{
  InputParameters params = DiracKernel::validParams();
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

AnisotropicNeutronPointSource::AnisotropicNeutronPointSource(const InputParameters & parameters)
  : DiracKernel(parameters)
  , _directions(getADMaterialProperty<std::vector<RealVectorValue>>("directions"))
  , _ordinate_index(getParam<unsigned int>("ordinate_index"))
  , _source_intensity(getParam<Real>("intensities"))
  , _angular_distribution(getFunction("phase_function"))
  , _source_location(getParam<Point>("points"))
{ }

void
AnisotropicNeutronPointSource::addPoints()
{
  addPoint(_source_location);
}

Real
AnisotropicNeutronPointSource::computeQpResidual()
{
  Point temp(MetaPhysicL::raw_value(_directions[_qp][_ordinate_index](0)),
             MetaPhysicL::raw_value(_directions[_qp][_ordinate_index](1)),
             MetaPhysicL::raw_value(_directions[_qp][_ordinate_index](2)));

  // Hijacking the MOOSE function system so the user can parse in an analytical
  // phase function for this point source.
  return -1.0 * _test[_i][_qp] * _angular_distribution.value(_t, temp)
              * _source_intensity;
}
