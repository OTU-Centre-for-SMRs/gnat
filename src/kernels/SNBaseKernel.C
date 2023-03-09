#include "SNBaseKernel.h"

InputParameters
SNBaseKernel::validParams()
{
  auto params = Kernel::validParams();
  params.addClassDescription("Provides basic functionality for neutron "
                             "transport kernels that require angular "
                             "quadrature sets. This kernel does NOT implement "
                             "computeQpResidual().");
  params.addRequiredRangeCheckedParam<unsigned int>("n_l",
                                                    "n_l > 0",
                                                    "Order of the polar Gauss-"
                                                    "Legendre quadrature set.");
  params.addRequiredRangeCheckedParam<unsigned int>("n_c",
                                                    "n_c > 0",
                                                    "Order of the azimuthal "
                                                    "Gauss-Chebyshev "
                                                    "quadrature set.");
  params.addParam<MooseEnum>("major_axis",
                             MooseEnum("x y z", "x"),
                             "Major axis of the angular quadrature. Allows the "
                             "polar angular quadrature to align with a cartesian "
                             "axis with minimal heterogeneity. Default is the "
                             "x-axis. This parameter is ignored for 1D and 2D "
                             "problems.");
  params.addParam<std::string>(
      "transport_system",
      "",
      "Name of the transport system which will consume the provided material properties. If one is "
      "not provided the first transport system will be used.");
  params.addRequiredParam<MooseEnum>("dimensionality",
                                     MooseEnum("1D_cartesian 2D_cartesian 3D_cartesian"),
                                     "Dimensionality and the coordinate system of the "
                                     "problem.");

  return params;
}

SNBaseKernel::SNBaseKernel(const InputParameters & parameters)
  : Kernel(parameters),
    _quadrature_set(getParam<unsigned int>("n_c"),
                    getParam<unsigned int>("n_l"),
                    getParam<MooseEnum>("major_axis").getEnum<MajorAxis>(),
                    getParam<MooseEnum>("dimensionality").getEnum<ProblemType>()),
    _symmetry_factor(1.0)
{
  switch (_quadrature_set.getProblemType())
  {
    case ProblemType::Cartesian1D:
      _symmetry_factor = 2.0 * M_PI;
      break;

    case ProblemType::Cartesian2D:
      _symmetry_factor = 2.0;
      break;

    case ProblemType::Cartesian3D:
      _symmetry_factor = 1.0;
      break;

    default:
      _symmetry_factor = 1.0;
      break;
  }
}

void
SNBaseKernel::cartesianToSpherical(const RealVectorValue & ordinate, Real & mu, Real & omega)
{
  switch (_quadrature_set.getAxis())
  {
    case MajorAxis::X:
      mu = ordinate(0);
      omega = std::acos(ordinate(1) / std::sqrt(1.0 - (mu * mu)));

      break;

    case MajorAxis::Y:
      mu = ordinate(1);
      omega = std::acos(ordinate(2) / std::sqrt(1.0 - (mu * mu)));

      break;

    case MajorAxis::Z:
      mu = ordinate(2);
      omega = std::acos(ordinate(0) / std::sqrt(1.0 - (mu * mu)));

      break;
  }
}
