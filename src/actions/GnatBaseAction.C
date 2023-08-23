#include "GnatBaseAction.h"

#include "FEProblemBase.h"
#include "AddVariableAction.h"
#include "InputParameterWarehouse.h"

InputParameters
GnatBaseAction::validParams()
{
  auto params = Action::validParams();
  params.addClassDescription("This action serves as a base for all of Gnat's "
                             "actions, and includes functionality and "
                             "parameters required to setup coupled transport + dispersion studies. "
                             "This action does not implement act().");
  params.addParam<MooseEnum>("debug_verbosity",
                             MooseEnum("level0 level1", "level1"),
                             "How verbose the debug output of the transport "
                             "system should be. level0 is fully verbose. "
                             "level1 outputs less debugging information.");
  params.addParam<std::vector<SubdomainName>>("block",
                                              "The list of blocks (ids or "
                                              "names) that this variable will "
                                              "be applied.");

  return params;
}

GnatBaseAction::GnatBaseAction(const InputParameters & params)
  : Action(params),
    _exec_type(ExecutionType::SteadySource),
    _debug_level(getParam<MooseEnum>("debug_verbosity").getEnum<DebugVerbosity>()),
    _p_type(ProblemType::Cartesian1D)
{
}

void
GnatBaseAction::initializeBase()
{
  // Grab all the subdomain IDs that the derived actions should be
  // applied to.
  for (const auto & block_name : getParam<std::vector<SubdomainName>>("block"))
    _subdomain_ids.insert(_problem->mesh().getSubdomainID(block_name));

  // Find out what sort of problem we're working with. Error if it's not a
  // cartesian coordinate system.
  for (const SubdomainID & id : _subdomain_ids)
  {
    if (_problem->getCoordSystem(id) != Moose::COORD_XYZ)
      mooseError("Transport simulations currently do not support "
                 "non-cartesian coordinate systems.");
  }

  // Set the execution type.
  if (_problem->isTransient())
    _exec_type = ExecutionType::Transient;

  switch (_mesh->dimension())
  {
    case 1u:
      _p_type = ProblemType::Cartesian1D;
      break;

    case 2u:
      _p_type = ProblemType::Cartesian2D;
      break;

    case 3u:
      _p_type = ProblemType::Cartesian3D;
      break;

    default:
      mooseError("Unknown mesh dimensionality.");
      break;
  }
}

void
GnatBaseAction::debugOutput(const std::string & level0,
                            const std::string & level1,
                            const std::string & color)
{
  switch (_debug_level)
  {
    case DebugVerbosity::Level0:
      if (level0 != "")
      {
        if (level0 == "\n")
          _console << COLOR_DEFAULT << std::endl;
        else
          _console << color << level0 << COLOR_DEFAULT << std::endl;
      }
      break;

    case DebugVerbosity::Level1:
      if (level1 != "")
      {
        if (level1 == "\n")
          _console << COLOR_DEFAULT << std::endl;
        else
          _console << color << level1 << COLOR_DEFAULT << std::endl;
      }
      break;

    default:
      break;
  }
}

void
GnatBaseAction::addVariable(const std::string & var_name)
{
  auto fe_type = AddVariableAction::feType(_pars);
  auto type = AddVariableAction::variableType(fe_type, false, false);
  auto var_params = _factory.getValidParams(type);
  var_params.applySpecificParameters(_pars, {"family", "order"});
  var_params.set<std::vector<Real>>("scaling") = {getParam<Real>("scaling")};

  if (_subdomain_ids.empty())
    _problem->addVariable(type, var_name, var_params);
  else
  {
    for (const SubdomainID & id : _subdomain_ids)
      var_params.set<std::vector<SubdomainName>>("block").push_back(Moose::stringify(id));

    _problem->addVariable(type, var_name, var_params);
  }

  debugOutput("      - Adding variable " + var_name + ".");
}
