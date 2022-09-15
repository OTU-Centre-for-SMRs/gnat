#pragma once

#include "ADKernel.h"

class ADIsotopeBase : public ADKernel
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

  ADIsotopeBase(const InputParameters & parameters);

protected:
  constexpr static Real _avogadros_number = 6.02214076e23;

  // Helper function to fetch the velocity at _qp;
  ADRealVectorValue getQpVelocity();

  // Functions for computing SUPG stabilization parameters.
  ADReal computeQpTau();
  ADReal computeQpTests();

  // Functions for converting between number and mass densities.
  ADReal adComputeNumberDensity(ADReal mass_density);
  ADReal adComputeMassDensity(ADReal number_density);

  // Enums to make reading the source code easier.
  enum class VelocityType
  {
    Constant = 0u,
    Function = 1u,
    VariableComponent = 2u,
    VariableCombined = 3u
  } _vel_type;

  enum class DensityType
  {
    NumberDensity = 0u,
    MassDensity = 1u
  } _primal_variable_type;

  // Molar mass to facilitate conversion between mass density and number density.
  const Real _molar_mass;

  // Velocity vectors for constant, function and variable velocity fields.
  std::vector<const ADVariableValue *> _variable_comp_vel;
  std::vector<const ADVectorVariableValue *> _variable_vec_vel;
  std::vector<const Function *> _function_vel;
  const RealVectorValue _constant_vel;

  const unsigned int _mesh_dims;
}; // class ADIsotopeBase
