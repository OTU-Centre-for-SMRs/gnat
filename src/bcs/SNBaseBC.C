#include "SNBaseBC.h"

InputParameters
SNBaseBC::validParams()
{
  auto params = IntegratedBC::validParams();
  params.addClassDescription("Provides basic functionality for the neutron "
                             "transport boundary conditions. This BC does NOT "
                             "implement computeQpResidual().");
  params.addRequiredParam<UserObjectName>(
      "aq", "The name of the angular quadrature provider user object.");

  return params;
}

SNBaseBC::SNBaseBC(const InputParameters & parameters)
  : IntegratedBC(parameters), _aq(getUserObject<AQProvider>("aq")), _symmetry_factor(1.0)
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
SNBaseBC::cartesianToSpherical(const RealVectorValue & ordinate, Real & mu, Real & omega)
{
  mu = 0.0;
  omega = 0.0;

  switch (_aq.getAxis())
  {
    case MajorAxis::X:
      mu = ordinate(0);

      if (ordinate(1) > 0.0)
      {
        if (ordinate(2) > 0.0)
          omega += std::atan(std::abs(ordinate(2)) / std::abs(ordinate(1)));
        else if (ordinate(2) < 0.0)
          omega += std::atan(std::abs(ordinate(2)) / std::abs(ordinate(1))) + libMesh::pi;
      }
      else if (ordinate(1) < 0.0)
      {
        if (ordinate(2) > 0.0)
          omega +=
              std::atan(std::abs(ordinate(2)) / std::abs(ordinate(1))) + 3.0 * libMesh::pi / 2.0;
        else if (ordinate(2) < 0.0)
          omega += std::atan(std::abs(ordinate(2)) / std::abs(ordinate(1))) + libMesh::pi / 2.0;
        else
          omega += libMesh::pi;
      }
      else
      {
        if (ordinate(2) > 0.0)
          omega += libMesh::pi / 2.0;
        else if (ordinate(2) < 0.0)
          omega += 3.0 * libMesh::pi / 2.0;
      }
      break;

    default:
      mooseError("This logic has not been implemented yet.");
  }
}
