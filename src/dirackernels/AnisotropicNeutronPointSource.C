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
  params.addRequiredRangeCheckedParam<unsigned int>("n_l",
                                                    "n_l > 0",
                                                    "Order of the polar Gauss-"
                                                    "Legendre quadrature set.");
  params.addRequiredRangeCheckedParam<unsigned int>("n_c",
                                                    "n_c > 0",
                                                    "Order of the azimuthal "
                                                    "Gauss-Chebyshev "
                                                    "quadrature set.");
  params.addParam<MooseEnum>("major_axis", MooseEnum("x y z", "x"),
                             "Major axis of the angular quadrature. Allows the "
                             "polar angular quadrature to align with a cartesian "
                             "axis with minimal heterogeneity. Default is the "
                             "x-axis. This parameter is ignored for 1D and 2D "
                             "problems.");
  params.addRequiredParam<MooseEnum>("dimensionality",
                                     MooseEnum("1D_cartesian 2D_cartesian 3D_cartesian"),
                                     "Dimensionality and the coordinate system of the "
                                     "problem.");
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
  , _quadrature_set(getParam<unsigned int>("n_c"),
                    getParam<unsigned int>("n_l"),
                    getParam<MooseEnum>("major_axis").getEnum<MajorAxis>(),
                    getParam<MooseEnum>("dimensionality").getEnum<ProblemType>())
  , _ordinate_index(getParam<unsigned int>("ordinate_index"))
  , _source_intensity(getParam<Real>("intensities"))
  , _angular_distribution(getFunction("phase_function"))
  , _source_location(getParam<Point>("points"))
  , _symmetry_factor(1.0)
{
  switch (_quadrature_set.getProblemType())
  {
    case ProblemType::Cartesian1D: _symmetry_factor = 2.0 * M_PI; break;
    case ProblemType::Cartesian2D: _symmetry_factor = 2.0; break;
    case ProblemType::Cartesian3D: _symmetry_factor = 1.0; break;
    default: _symmetry_factor = 1.0; break;
  }
}

void
AnisotropicNeutronPointSource::addPoints()
{
  addPoint(_source_location);
}

Real
AnisotropicNeutronPointSource::computeQpResidual()
{
  Point temp(_quadrature_set.direction(_ordinate_index)(0),
             _quadrature_set.direction(_ordinate_index)(1),
             _quadrature_set.direction(_ordinate_index)(2));

  // Hijacking the MOOSE function system so the user can parse in an analytical
  // phase function for this point source.
  return (-1.0 / M_PI) * _test[_i][_qp] * _angular_distribution.value(_t, temp)
         * _source_intensity * _symmetry_factor;
}
