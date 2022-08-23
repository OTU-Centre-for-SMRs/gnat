#include "ADDGIsotopeConstAdvection.h"

registerMooseObject("GnatApp", ADDGIsotopeConstAdvection);

InputParameters
ADDGIsotopeConstAdvection::validParams()
{
  auto params = ADDGKernel::validParams();
  params.addClassDescription("Computes the discontinuous interior face term "
                             "for the isotope scalar transport equation with a "
                             "constant velocity vector using upwinding for the "
                             "numerical flux.");

  params.addRequiredParam<RealVectorValue>("velocity", "The velocity vector.");

  return params;
}

ADDGIsotopeConstAdvection::ADDGIsotopeConstAdvection(const InputParameters & parameters)
  : ADDGKernel(parameters)
  , _vel(getParam<RealVectorValue>("velocity"))
{ }

ADReal
ADDGIsotopeConstAdvection::computeQpResidual(Moose::DGResidualType type)
{

  ADReal res = 0.0;
  ADReal n_dot_v = _vel * _normals[_qp];
  switch (type)
  {
    case Moose::Element:
      if (n_dot_v >= 0)
        res += n_dot_v * _u[_qp] * _test[_i][_qp];
      else
        res += n_dot_v * _u_neighbor[_qp] * _test[_i][_qp];
      break;

    case Moose::Neighbor:
      if (n_dot_v >= 0)
        res -= n_dot_v * _u[_qp] * _test_neighbor[_i][_qp];
      else
        res -= n_dot_v * _u_neighbor[_qp] * _test_neighbor[_i][_qp];
      break;
  }

  return res;
}
