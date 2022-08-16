#include "CommonGnatAction.h"

#include "ActionWarehouse.h"

#include "GnatBaseAction.h"

registerMooseAction("GnatApp", CommonGnatAction, "meta_action");

InputParameters
CommonGnatAction::validParams()
{
  auto params = GnatBaseAction::validParams();
  params.addClassDescription("Action to store common Gnat parameters.");

  return params;
}

CommonGnatAction::CommonGnatAction(const InputParameters & params)
  : Action(params)
{ }

void
CommonGnatAction::act()
{
  // Check if sub-blocks block are found which will use the common parameters.
  auto action = _awh.getActions<GnatBaseAction>();
  if (action.size() == 0)
    mooseWarning("Common parameters are supplied, but not used in ",
                 parameters().blockLocation(), ".");
}
