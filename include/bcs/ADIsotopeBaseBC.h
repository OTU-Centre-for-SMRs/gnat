#pragma once

#include "ADIntegratedBC.h"

class ADIsotopeBaseBC : public ADIntegratedBC
{
public:
  // An enum to get around MooseEnum type checking for enums.
  enum class MooseEnumVelocityType
  {
    Constant = 0u,
    Function = 1u,
    Variable = 2u
  };

  static InputParameters validParams();

  ADIsotopeBaseBC(const InputParameters & parameters);

protected:
  // Helper function to fetch the velocity at _qp;
  ADRealVectorValue getQpVelocity();

  // Functions for computing SUPG stabilization parameters.
  ADReal computeQpTau();
  ADReal computeQpTests();

  // Enums to make reading the source code easier.
  enum class VelocityType
  {
    Constant = 0u,
    Function = 1u,
    VariableComponent = 2u,
    VariableCombined = 3u
  } _vel_type;

  // Velocity vectors for constant, function and variable velocity fields.
  std::vector<const ADVariableValue *> _variable_comp_vel;
  std::vector<const ADVectorVariableValue *> _variable_vec_vel;
  std::vector<const Function *> _function_vel;
  const RealVectorValue _constant_vel;

  const unsigned int _mesh_dims;
};
