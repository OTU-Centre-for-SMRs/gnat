#include "ADIsotopeVarAdvection.h"

#include "MooseMesh.h"

registerMooseObject("GnatApp", ADIsotopeVarAdvection);

InputParameters
ADIsotopeVarAdvection::validParams()
{
  auto params = ADKernel::validParams();
  params.addClassDescription("Computes the advection term for the isotope "
                             "scalar transport equation with a coupled "
                             "velocity field: "
                             "$(- ( \\vec{\\nabla}\\psi_{j}\\cdot\\vec{v},"
                             "N_{i} )_{V}$. This kernel is unstable and can "
                             "only be used with discontinuous finite "
                             "elements.");
  params.addCoupledVar("u", "The x-component of the velocity field.");
  params.addCoupledVar("v", "The v-component of the velocity field.");
  params.addCoupledVar("w", "The w-component of the velocity field.");

  params.addCoupledVar("velocity",
                       "A vector velocity value as opposed to using "
                       "individual velocity components.");

  return params;
}

ADIsotopeVarAdvection::ADIsotopeVarAdvection(const InputParameters & parameters)
  : ADKernel(parameters)
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
ADIsotopeVarAdvection::computeQpResidual()
{
  ADReal res = 0.0;

  if (!isCoupled("velocity"))
  {
    res -= _mesh_dims >= 1u ? _grad_test[_i][_qp](0) * _u_vel[_qp] : 0.0;
    res -= _mesh_dims >= 2u ? _grad_test[_i][_qp](1) * _v_vel[_qp] : 0.0;
    res -= _mesh_dims == 3u ? _grad_test[_i][_qp](2) * _w_vel[_qp] : 0.0;
    res *= _u[_qp];
  }
  else
    res -= _grad_test[_i][_qp] * _vel[_qp] * _u[_qp];

  return res;
}
