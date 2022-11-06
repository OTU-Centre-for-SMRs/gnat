#include "SetupIsotopeSystemAction.h"

#include "ActionWarehouse.h"

#include "AddMobileIsotopeAction.h"

registerMooseAction("GnatApp", SetupIsotopeSystemAction, "meta_action");

InputParameters
SetupIsotopeSystemAction::validParams()
{
  auto params = Action::validParams();
  params.addClassDescription("Action to store common isotope mass transport parameters.");
  //----------------------------------------------------------------------------
  // Parameters required for mass transport and stabilization. These are applied
  // to all isotope kernels and BCs.
  params.addRequiredParam<std::vector<VariableName>>(
      "isotopes", "A list of all isotopes considered in the isotope system, mobile or not.");
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

SetupIsotopeSystemAction::SetupIsotopeSystemAction(const InputParameters & parameters)
  : Action(parameters)
{
}

void
SetupIsotopeSystemAction::act()
{
  // Check if sub-blocks block are found which will use the common parameters.
  auto action = _awh.getActions<AddMobileIsotopeAction>();
  if (action.size() == 0)
    mooseWarning("Velocity field parameters are supplied, but not used in ",
                 parameters().blockLocation(),
                 ".");
}
