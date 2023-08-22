#include "SNBaseDiracKernel.h"

InputParameters
SNBaseDiracKernel::validParams()
{
  auto params = DiracKernel::validParams();
  params.addClassDescription("Provides basic functionality for neutron "
                             "transport Dirac kernels that require angular "
                             "quadrature sets. This kernel does NOT implement "
                             "computeQpResidual() or computeQpJacobian().");
  params.addRequiredParam<UserObjectName>(
      "aq", "The name of the angular quadrature provider user object.");

  return params;
}

SNBaseDiracKernel::SNBaseDiracKernel(const InputParameters & parameters)
  : DiracKernel(parameters), _aq(getUserObject<AQProvider>("aq")), _symmetry_factor(1.0)
{
  switch (_aq.getProblemType())
  {
    case ProblemType::Cartesian1D:
      _symmetry_factor = 2.0 * libMesh::pi;
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
SNBaseDiracKernel::cartesianToSpherical(const RealVectorValue & ordinate, Real & mu, Real & omega)
{
  switch (_aq.getAxis())
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
