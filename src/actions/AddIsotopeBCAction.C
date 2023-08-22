#include "AddIsotopeBCAction.h"

#include "FEProblem.h"
#include "BoundaryCondition.h"
#include "InputParameterWarehouse.h"

#include "ADIsotopeBaseBC.h"
#include "ADIsotopeBase.h"

#include "MobileDepletionSystemAction.h"

registerMooseAction("GnatApp", AddIsotopeBCAction, "add_bc");

InputParameters
AddIsotopeBCAction::validParams()
{
  auto params = MooseObjectAction::validParams();
  // params += GnatBaseAction::validParams();

  params.addClassDescription("Adds a boundary condition to all isotopes in the isotope system.");
  params.addParam<std::vector<VariableName>>(
      "excluded_nuclides",
      std::vector<VariableName>(),
      "Isotopes in the master list that this boundary condition should "
      "not be applied to.");

  return params;
}

AddIsotopeBCAction::AddIsotopeBCAction(const InputParameters & parameters)
  : MooseObjectAction(parameters),
    _exclude(getParam<std::vector<VariableName>>("excluded_nuclides"))
{
}

void
AddIsotopeBCAction::act()
{
  const unsigned int mesh_dims = _problem->mesh().dimension();

  const auto nuclide_actions = _awh.getActions<MobileDepletionSystemAction>();
  if (nuclide_actions.size() != 1u)
    mooseError("There can only be a single MobileDepletionSystemAction at a time. There are currently " +
               Moose::stringify(nuclide_actions.size() + " MobileDepletionSystemActions."));
  const auto & nuclide_action = (*nuclide_actions[0u]);

  if (mesh_dims >= 2u && !nuclide_action.isParamValid("v"))
    mooseError("In 2D or 3D the v component of the velocity must be supplied using the 'v' "
               "parameter.");
  if (mesh_dims >= 3u && !nuclide_action.isParamValid("w"))
    mooseError("In 3D, the w component of the velocity must be supplied using the 'w' parameter.");

  for (const auto & var : nuclide_action.getParam<std::vector<VariableName>>("nuclides"))
  {
    auto r = std::find(_exclude.begin(), _exclude.end(), var);
    if (r != _exclude.end())
      continue;

    _moose_object_pars.set<NonlinearVariableName>("variable") = var;

    _moose_object_pars.set<MooseFunctorName>("u") = nuclide_action.getParam<MooseFunctorName>("u");
    if (nuclide_action.isParamValid("v"))
      _moose_object_pars.set<MooseFunctorName>("v") =
          nuclide_action.getParam<MooseFunctorName>("v");
    if (nuclide_action.isParamValid("w"))
      _moose_object_pars.set<MooseFunctorName>("v") =
          nuclide_action.getParam<MooseFunctorName>("w");

    _problem->addBoundaryCondition(_type, _name + "_" + var, _moose_object_pars);
  }
}
