#include "AddIsotopeBCAction.h"

#include "FEProblem.h"
#include "BoundaryCondition.h"
#include "InputParameterWarehouse.h"

#include "SetupNuclideSystemAction.h"
#include "ADIsotopeBaseBC.h"
#include "ADIsotopeBase.h"

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
      "not be applied to..");

  // Add common isotope system parameters.
  params += SetupNuclideSystemAction::validParams();
  params.makeParamNotRequired<MooseEnum>("velocity_type");
  params.makeParamNotRequired<std::vector<VariableName>>("nuclides");

  return params;
}

AddIsotopeBCAction::AddIsotopeBCAction(const InputParameters & parameters)
  : MooseObjectAction(parameters),
    _p_type(ProblemType::Cartesian1D),
    _master_isotope_list(getParam<std::vector<VariableName>>("nuclides")),
    _exclude(getParam<std::vector<VariableName>>("excluded_nuclides"))
{
  // Check if a container block exists with isotope parameters. If yes, apply them.
  // FIX THIS: The most janky way to fix this breaking MOOSE change.
  auto isotope_system_actions = _awh.getActions<SetupNuclideSystemAction>();
  if (isotope_system_actions.size() == 1)
  {
    const auto & params = _app.getInputParameterWarehouse().getInputParameters();
    InputParameters & pars(*(params.find(uniqueActionName())->second.get()));
    pars.applyParameters(isotope_system_actions[0]->parameters());
  }
}

void
AddIsotopeBCAction::applyIsotopeParameters(InputParameters & params)
{
  params.set<MooseEnum>("velocity_type") = getParam<MooseEnum>("velocity_type");

  switch (getParam<MooseEnum>("velocity_type").getEnum<ADIsotopeBase::MooseEnumVelocityType>())
  {
    case ADIsotopeBase::MooseEnumVelocityType::Constant:
      params.set<RealVectorValue>("constant_velocity") =
          getParam<RealVectorValue>("constant_velocity");
      break;

    case ADIsotopeBase::MooseEnumVelocityType::Function:
      switch (_p_type)
      {
        case ProblemType::Cartesian1D:
          params.set<FunctionName>("u_function") = getParam<FunctionName>("u_function");
          break;

        case ProblemType::Cartesian2D:
          params.set<FunctionName>("u_function") = getParam<FunctionName>("u_function");
          params.set<FunctionName>("v_function") = getParam<FunctionName>("v_function");
          break;

        case ProblemType::Cartesian3D:
          params.set<FunctionName>("u_function") = getParam<FunctionName>("u_function");
          params.set<FunctionName>("v_function") = getParam<FunctionName>("v_function");
          params.set<FunctionName>("w_function") = getParam<FunctionName>("w_function");
          break;

        default:
          break;
      }
      break;

    case ADIsotopeBase::MooseEnumVelocityType::Variable:
      if (isParamValid("u_var"))
      {
        switch (_p_type)
        {
          case ProblemType::Cartesian1D:
            params.set<std::vector<VariableName>>("u_var").emplace_back(
                getParam<VariableName>("u_var"));

            break;

          case ProblemType::Cartesian2D:
            params.set<std::vector<VariableName>>("u_var").emplace_back(
                getParam<VariableName>("u_var"));
            params.set<std::vector<VariableName>>("v_var").emplace_back(
                getParam<VariableName>("v_var"));

            break;

          case ProblemType::Cartesian3D:
            params.set<std::vector<VariableName>>("u_var").emplace_back(
                getParam<VariableName>("u_var"));
            params.set<std::vector<VariableName>>("v_var").emplace_back(
                getParam<VariableName>("v_var"));
            params.set<std::vector<VariableName>>("w_var").emplace_back(
                getParam<VariableName>("w_var"));

            break;

          default:
            break;
        }
      }
      else
        params.set<std::vector<VariableName>>("vector_velocity")
            .emplace_back(getParam<VariableName>("vector_velocity"));
      break;
  }
}

void
AddIsotopeBCAction::act()
{
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

  for (const auto & var : _master_isotope_list)
  {
    {
      auto r = std::find(_exclude.begin(), _exclude.end(), var);
      if (r != _exclude.end())
        continue;
    }

    _moose_object_pars.set<NonlinearVariableName>("variable") = var;
    applyIsotopeParameters(_moose_object_pars);

    _problem->addBoundaryCondition(_type, _name + "_" + var, _moose_object_pars);
  }
}
