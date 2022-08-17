#include "GnatBaseAction.h"

#include "FEProblemBase.h"

#include "CommonGnatAction.h"

InputParameters
GnatBaseAction::validParams()
{
  auto params = Action::validParams();
  params.addClassDescription("This action serves as a base for all of Gnat's "
                             "actions, and includes functionality and "
                             "parameters required to setup neutron activation "
                             "studies. This action does not implement act().");
  params.addParam<std::string>("flux_moment_names",
                               "flux_moment",
                               "Variable names for the moments of the angular "
                               "flux. The output format for the group flux "
                               "moments will be of the form "
                               "{flux_moment_names}_g_l_m.");
  params.addParam<unsigned int>("num_groups",
                                "The number of spectral energy groups in the "
                                "problem.");
  params.addParam<unsigned int>("max_anisotropy", 0,
                                "The maximum degree of anisotropy to evaluate. "
                                "Defaults to 0 for isotropic scattering.");
  params.addParam<MooseEnum>("execution_type",
                             MooseEnum("steady transient"),
                             "The method of execution for the problem. Options "
                             "are steady-state source driven problems and "
                             "transient source problems.");
  params.addParam<MooseEnum>("debug_verbosity",
                             MooseEnum("level0 level1", "level1"),
                             "How verbose the debug output of the transport "
                             "system should be. level0 is fully verbose. "
                             "level1 outputs less debugging information.");
  params.addParam<std::vector<SubdomainName>>("block",
                                              "The list of blocks (ids or "
                                              "names) that this variable will "
                                              "be applied.");

  params.addParamNamesToGroup("num_groups execution_type", "Required");

  return params;
}

GnatBaseAction::GnatBaseAction(const InputParameters & params)
  : Action(params)
  , _num_groups(0u)
  , _max_eval_anisotropy(0u)
  , _exec_type(ExecutionType::SteadySource)
  , _flux_moment_name(getParam<std::string>("flux_moment_names"))
  , _debug_level(getParam<MooseEnum>("debug_verbosity").getEnum<DebugVerbosity>())
{
  // Check if a container block exists with common parameters. If yes, apply them.
  auto common_actions = _awh.getActions<CommonGnatAction>();
  if (common_actions.size() == 1)
    _pars.applyParameters(common_actions[0]->parameters());

  // Error check to make sure the user has supplied mandatory parameters.
  if (!isParamValid("num_groups"))
  {
    mooseError("The number of neutron energy groups has not been supplied. "
               "Please initialize num_groups to a value greater than 0.");
  }
  if (getParam<unsigned int>("num_groups") == 0u)
  {
    mooseError("Invalid number of neutron energy groups has been supplied. "
               "num_groups must be greater than 0. ");
  }
  if (!isParamValid("execution_type"))
  {
    mooseError("The problem execution type must be provided. Please initialize "
               "execution_type to steady or transient.");
  }

  _num_groups = getParam<unsigned int>("num_groups");
  _max_eval_anisotropy = getParam<unsigned int>("max_anisotropy");
  _exec_type = getParam<MooseEnum>("execution_type").getEnum<ExecutionType>();
  _flux_moment_name = getParam<std::string>("flux_moment_names");
  _debug_level = getParam<MooseEnum>("debug_verbosity").getEnum<DebugVerbosity>();
}

void
GnatBaseAction::initializeBase()
{
  // Grab all the subdomain IDs that the neutron transport action should be
  // applied to.
  for (const auto & block_name : getParam<std::vector<SubdomainName>>("block"))
    _subdomain_ids.insert(_problem->mesh().getSubdomainID(block_name));

  // Find out what sort of problem we're working with. Error if it's not a
  // cartesian coordinate system.
  for (const SubdomainID & id : _subdomain_ids)
  {
    if (_problem->getCoordSystem(id) != Moose::COORD_XYZ)
      mooseError("Neutron transport simulations currently do not support "
                 "non-cartesian coordinate systems.");
  }

  // Set the enum so quadrature sets can be determined appropriately.
  unsigned int _num_group_moments = 0u;
  switch (_mesh->dimension())
  {
    case 1u:
      _p_type = ProblemType::Cartesian1D;
      _num_group_moments = (_max_eval_anisotropy + 1u);

      for (unsigned int g = 0; g < _num_groups; ++g)
      {
        // Set up variable names for the group flux moments.
        _group_flux_moments.emplace(g, std::vector<VariableName>());
        _group_flux_moments[g].reserve(_num_group_moments);
        for (unsigned int l = 0; l <= _max_eval_anisotropy; ++l)
        {
          _group_flux_moments[g].emplace_back(_flux_moment_name + "_"
                                              + Moose::stringify(g + 1u) + "_"
                                              + Moose::stringify(l) + "_"
                                              + Moose::stringify(0));
        }
      }
      break;

    case 2u:
      _p_type = ProblemType::Cartesian2D;
      _num_group_moments = (_max_eval_anisotropy + 1u) * (_max_eval_anisotropy + 2u) / 2u;

      for (unsigned int g = 0; g < _num_groups; ++g)
      {
        // Set up variable names for the group flux moments.
        _group_flux_moments.emplace(g, std::vector<VariableName>());
        _group_flux_moments[g].reserve(_num_group_moments);
        for (unsigned int l = 0; l <= _max_eval_anisotropy; ++l)
        {
          for (int m = 0; m <= static_cast<int>(l); ++m)
          {
            _group_flux_moments[g].emplace_back(_flux_moment_name + "_"
                                                + Moose::stringify(g + 1u) + "_"
                                                + Moose::stringify(l) + "_"
                                                + Moose::stringify(m));
          }
        }
      }
      break;

    case 3u:
      _p_type = ProblemType::Cartesian3D;
      _num_group_moments = (_max_eval_anisotropy + 1u) * (_max_eval_anisotropy + 1u);

      for (unsigned int g = 0; g < _num_groups; ++g)
      {
        // Set up variable names for the group flux moments.
        _group_flux_moments.emplace(g, std::vector<VariableName>());
        _group_flux_moments[g].reserve(_num_group_moments);
        for (unsigned int l = 0; l <= _max_eval_anisotropy; ++l)
        {
          for (int m = -1 * static_cast<int>(l); m <= static_cast<int>(l); ++m)
          {
            _group_flux_moments[g].emplace_back(_flux_moment_name + "_"
                                                + Moose::stringify(g + 1u) + "_"
                                                + Moose::stringify(l) + "_"
                                                + Moose::stringify(m));
          }
        }
      }
      break;

    default:
      mooseError("Unknown mesh dimensionality.");
      break;
  }
}

void
GnatBaseAction::debugOutput(const std::string & level0, const std::string & level1)
{
  switch (_debug_level)
  {
    case DebugVerbosity::Level0:
      if (level0 != "")
      {
        if (level0 == "\n")
          _console << std::endl;
        else
          _console << level0 << std::endl;
      }
      break;

    case DebugVerbosity::Level1:
      if (level1 != "")
      {
        if (level1 == "\n")
          _console << std::endl;
        else
          _console << level1 << std::endl;
      }
      break;

    default:
      break;
  }
}
