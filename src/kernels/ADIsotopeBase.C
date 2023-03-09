#include "ADIsotopeBase.h"

#include "Function.h"

InputParameters
ADIsotopeBase::validParams()
{
  auto params = ADKernel::validParams();
  params.addClassDescription("A base class whick computes the S/U-PG "
                             "stabilization term for the isotope mass "
                             "transport equation.");
  params.addRequiredParam<MooseEnum>("velocity_type",
                                     MooseEnum("constant function variable"),
                                     "An indicator for which type of velocity "
                                     "field should be used.");
  params.addParam<RealVectorValue>(
      "constant_velocity", RealVectorValue(0.0), "A constant velocity field.");
  params.addParam<FunctionName>("u_function",
                                "The x-component of the function "
                                "velocity field.");
  params.addParam<FunctionName>("v_function",
                                "The y-component of the function "
                                "velocity field.");
  params.addParam<FunctionName>("w_function",
                                "The z-component of the function "
                                "velocity field.");
  params.addCoupledVar("u_var",
                       "The x-component of the variable velocity "
                       "field.");
  params.addCoupledVar("v_var",
                       "The y-component of the variable velocity "
                       "field.");
  params.addCoupledVar("w_var",
                       "The z-component of the variable velocity "
                       "field.");
  params.addCoupledVar("vector_velocity",
                       "A vector variable velocity field as opposed to using "
                       "individual velocity components.");

  params.addCoupledVar("eddy_diffusivity",
                       "The eddy diffusivity computed with turbulence modelling.");
  // We need some ghosting for the finite volume fields (we use neighbor information to compute
  // gradients).
  params.addParam<unsigned short>("ghost_layers", 2, "The number of layers of elements to ghost.");
  params.addRelationshipManager(
      "ElementSideNeighborLayers",
      Moose::RelationshipManagerType::GEOMETRIC | Moose::RelationshipManagerType::ALGEBRAIC,
      [](const InputParameters & obj_params, InputParameters & rm_params)
      {
        rm_params.set<unsigned short>("layers") = obj_params.get<unsigned short>("ghost_layers");
        rm_params.set<bool>("use_displaced_mesh") = obj_params.get<bool>("use_displaced_mesh");
      });

  return params;
}

ADIsotopeBase::ADIsotopeBase(const InputParameters & parameters)
  : ADKernel(parameters),
    _vel_type(VelocityType::Constant),
    _constant_vel(getParam<RealVectorValue>("constant_velocity")),
    _mat_diff(getADMaterialProperty<Real>(
        "isotope_diff_" + Moose::stringify(getParam<NonlinearVariableName>("variable")))),
    _eddy_diffusivity(isParamValid("eddy_diffusivity") ? &coupledValue("eddy_diffusivity")
                                                       : nullptr),
    _mesh_dims(_fe_problem.mesh().dimension())
{
  switch (getParam<MooseEnum>("velocity_type").getEnum<MooseEnumVelocityType>())
  {
    case MooseEnumVelocityType::Constant:
      // Constant velocity.
      _vel_type = VelocityType::Constant;
      break;

    case MooseEnumVelocityType::Function:
      // Function velocity.
      _vel_type = VelocityType::Function;

      if (isParamValid("u_function"))
      {
        _function_vel.emplace_back(&getFunction("u_function"));

        if (isParamValid("v_function"))
        {
          _function_vel.emplace_back(&getFunction("v_function"));

          if (isParamValid("w_function"))
            _function_vel.emplace_back(&getFunction("w_function"));
        }
      }

      if (_function_vel.size() != _mesh_dims)
      {
        mooseError("The number of provided velocity component functions does "
                   "not match the mesh dimensionality.");
      }

      break;

    case MooseEnumVelocityType::Variable:
      // Vector variable velocity.
      if (isCoupled("vector_velocity"))
      {
        // Vector variable velocity.
        _vel_type = VelocityType::VariableCombined;

        _variable_vec_vel.emplace_back(&adCoupledVectorValue("vector_velocity"));
      }
      else
      {
        // Component variable velocity.
        _vel_type = VelocityType::VariableComponent;

        if (isCoupled("u_var"))
        {
          _variable_comp_vel.emplace_back(&adCoupledValue("u_var"));

          if (isCoupled("v_var"))
          {
            _variable_comp_vel.emplace_back(&adCoupledValue("v_var"));

            if (isCoupled("w_var"))
              _variable_comp_vel.emplace_back(&adCoupledValue("w_var"));
          }
        }

        if (_variable_comp_vel.size() != _mesh_dims)
        {
          mooseError("The number of provided velocity component variables does "
                     "not match the mesh dimensionality. Components: " +
                     Moose::stringify(_variable_comp_vel.size()) +
                     ". Mesh dimensions:  " + Moose::stringify(_mesh_dims) + ".");
        }
      }
      break;
  }
}

ADRealVectorValue
ADIsotopeBase::getQpVelocity()
{
  ADRealVectorValue qp_velocity = ADRealVectorValue(0.0, 0.0, 0.0);
  switch (_vel_type)
  {
    case VelocityType::Constant:
      qp_velocity = _constant_vel;
      break;

    case VelocityType::Function:
      for (unsigned int i = 0u; i < _mesh_dims; ++i)
        qp_velocity(i) = (*_function_vel[i]).value(_t, _q_point[_qp]);
      break;

    case VelocityType::VariableComponent:
      for (unsigned int i = 0u; i < _mesh_dims; ++i)
        qp_velocity(i) = (*_variable_comp_vel[i])[_qp];
      break;

    case VelocityType::VariableCombined:
      qp_velocity = (*_variable_vec_vel[0])[_qp];
      break;
  }

  return qp_velocity;
}

ADReal
ADIsotopeBase::computeQpTau()
{
  ADReal inv_tau_sq = 0.0;

  const auto h_min = _current_elem->hmin();
  const auto h2 = h_min * h_min;

  // Advection contribution.
  const auto speed =
      (getQpVelocity() + ADRealVectorValue(libMesh::TOLERANCE * libMesh::TOLERANCE)).norm();
  inv_tau_sq += 4.0 * speed * speed / h2;

  // Transient contribution.
  if (_is_transient)
    inv_tau_sq += 4.0 / (_dt * _dt);

  // Diffusive contribution.
  if (_eddy_diffusivity)
    inv_tau_sq += 4.0 * (*_eddy_diffusivity)[_qp] / h2;
  inv_tau_sq += 4.0 * _mat_diff[_qp] / h2;

  return 1.0 / std::sqrt(inv_tau_sq);
}

ADReal
ADIsotopeBase::computeQpTests()
{
  ADReal tau = computeQpTau();
  return _test[_i][_qp] + tau * getQpVelocity() * _grad_test[_i][_qp];
}
