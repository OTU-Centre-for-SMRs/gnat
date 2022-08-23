#include "ADDGIsotopeVarAdvection.h"

registerMooseObject("GnatApp", ADDGIsotopeVarAdvection);

InputParameters
ADDGIsotopeVarAdvection::validParams()
{
  auto params = ADDGKernel::validParams();
  params.addClassDescription("Computes the discontinuous interior face term "
                             "for the isotope scalar transport equation with a "
                             "coupled velocity field using upwinding for the "
                             "numerical flux.");
  params.addCoupledVar("u", "The x-component of the velocity field.");
  params.addCoupledVar("v", "The v-component of the velocity field.");
  params.addCoupledVar("w", "The w-component of the velocity field.");

  params.addCoupledVar("velocity",
                       "A vector velocity value as opposed to using "
                       "individual velocity components.");

  return params;
}

ADDGIsotopeVarAdvection::ADDGIsotopeVarAdvection(const InputParameters & parameters)
  : ADDGKernel(parameters)
  , _vel(isCoupled("velocity") ? adCoupledVectorValue("velocity") : _vel_temp)
  , _u_vel(isCoupled("u") ? adCoupledValue("u") : _u_vel_temp)
  , _v_vel(isCoupled("v") ? adCoupledValue("v") : _v_vel_temp)
  , _w_vel(isCoupled("w") ? adCoupledValue("w") : _w_vel_temp)
  , _mesh_dims(_fe_problem.mesh().dimension())
{
  if (_mesh_dims >= 1u && !isCoupled("u") && !isCoupled("velocity"))
    mooseError("Either velocity or u must be declared for all problems.");

  if (_mesh_dims >= 2u && !isCoupled("v") && !isCoupled("velocity"))
    mooseError("Either velocity or v must be declared for 2D and 3D problems.");

  if (_mesh_dims == 3u && !isCoupled("w") && !isCoupled("velocity"))
    mooseError("Either velocity or w must be declared for 3D problems.");
}

ADReal
ADDGIsotopeVarAdvection::computeQpResidual(Moose::DGResidualType type)
{
  ADRealVectorValue velocity = ADRealVectorValue(0.0);
  if (!isCoupled("velocity"))
  {
    velocity(0) = _u_vel[_qp];
    velocity(1) = _v_vel[_qp];
    velocity(2) = _w_vel[_qp];
  }
  else
    velocity = _vel[_qp];

  ADReal res = 0.0;
  ADReal n_dot_v = velocity * _normals[_qp];
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
