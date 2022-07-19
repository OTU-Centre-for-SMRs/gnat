#include "ADNeutronMatSourceBC.h"

registerMooseObject("GnatApp", ADNeutronMatSourceBC);

InputParameters
ADNeutronMatSourceBC::validParams()
{
  auto params = ADIntegratedBC::validParams();
  params.addClassDescription("Computes the surface source boundary condition "
                             "with the weak form given by "
                             "$\\langle \\psi_{j},\\, \\hat{n}\\cdot"
                             "\\hat{\\Omega}\\Psi_{inc,\\, g}\\rangle_{\\Gamma_{i}}$, "
                             "$\\hat{n}\\cdot\\hat{\\Omega} \\leq 0$. "
                             "This kernel should not be exposed to the user, "
                             "instead being enabled through a transport action.");
  params.addRequiredParam<MooseEnum>("dimensionality",
                                     "Dimensionality and the coordinate system of the "
                                     "problem.");
  params.addRequiredRangeCheckedParam<unsigned int>("ordinate_index",
                                                    "ordinate_index >= 0",
                                                    "The discrete ordinate index "
                                                    "of the current angular "
                                                    "flux.");

  return params;
}

ADNeutronMatSourceBC::ADNeutronMatSourceBC(const InputParameters & parameters)
  : ADIntegratedBC(parameters)
  , _type(getParam<MooseEnum>("dimensionality").getEnum<ProblemType>())
  , _ordinate_index(getParam<unsigned int>("ordinate_index"))
  , _directions(getMaterialProperty<std::vector<RealVectorValue>>("directions"))
  , _surface_source(getADMaterialProperty<std::vector<Real>>("surface_source"))
  , _symmetry_factor(1.0)
{
  switch (_type)
  {
    case ProblemType::Cartesian1D: _symmetry_factor = 2.0 * M_PI; break;
    case ProblemType::Cartesian2D: _symmetry_factor = 2.0; break;
    case ProblemType::Cartesian3D: _symmetry_factor = 1.0; break;
    default: _symmetry_factor = 1.0; break;
  }
}

ADReal
ADNeutronMatSourceBC::computeQpResidual()
{
  if (_ordinate_index >= _directions[_qp].size()
      || _ordinate_index >= _surface_source[_qp].size())
  {
    mooseError("The ordinates index exceeds the number of quadrature points "
               "and/or source values.");
  }

  ADReal res = 0.0;
  ADReal n_dot_omega = _directions[_qp][_ordinate_index] * _normals[_qp];
  if (n_dot_omega > 0.0)
    res += _test[_i][_qp] * _u[_qp] * n_dot_omega;
  else
  {
    res += _test[_i][_qp] * _surface_source[_qp][_ordinate_index] * n_dot_omega
           * _symmetry_factor;
  }

  return res;
}
