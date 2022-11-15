#include "SetupNuclideSystemAction.h"

#include "ActionWarehouse.h"

#include "AddMobileIsotopeAction.h"

registerMooseAction("GnatApp", SetupNuclideSystemAction, "meta_action");

InputParameters
SetupNuclideSystemAction::validParams()
{
  auto params = Action::validParams();
  params.addClassDescription("Action to store common isotope mass transport parameters.");
  //----------------------------------------------------------------------------
  // Parameters required for mass transport and stabilization. These are applied
  // to all nuclides kernels and BCs.
  params.addRequiredParam<std::vector<VariableName>>(
      "nuclides", "A list of all nuclides considered in the system, mobile or not.");
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

  return params;
}

SetupNuclideSystemAction::SetupNuclideSystemAction(const InputParameters & parameters)
  : Action(parameters)
{
}

void
SetupNuclideSystemAction::act()
{
}
