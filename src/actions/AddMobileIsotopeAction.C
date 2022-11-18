#include "AddMobileIsotopeAction.h"

#include "AddVariableAction.h"
#include "SetupNuclideSystemAction.h"
#include "InputParameterWarehouse.h"

#include "ADIsotopeBase.h"

registerMooseAction("GnatApp", AddMobileIsotopeAction, "add_variable");
registerMooseAction("GnatApp", AddMobileIsotopeAction, "add_ic");
registerMooseAction("GnatApp", AddMobileIsotopeAction, "add_kernel");
registerMooseAction("GnatApp", AddMobileIsotopeAction, "add_material");

InputParameters
AddMobileIsotopeAction::validParams()
{
  auto params = GnatBaseAction::validParams();
  params.addClassDescription("This action adds a variable which represents the "
                             "mass/number density of a mobile isotope, "
                             "alongside the kernels required to describe the "
                             "time evolution of the isotope.");

  //----------------------------------------------------------------------------
  // Parameters for variables.
  params.addRequiredParam<MooseEnum>("family",
                                     AddVariableAction::getNonlinearVariableFamilies(),
                                     "Specifies the family of FE shape functions to "
                                     "use for this variable.");
  params.addRequiredParam<MooseEnum>("order",
                                     AddVariableAction::getNonlinearVariableOrders(),
                                     "Specifies the order of the FE shape "
                                     "function to use for this variable "
                                     "(additional orders not listed are "
                                     "allowed).");
  params.addParam<Real>("scaling",
                        1.0,
                        "Specifies a scaling factor to apply to "
                        "this variable.");
  params.addRequiredParam<VariableName>("isotope_name", "The identifying name of the isotope.");

  //----------------------------------------------------------------------------
  // Properties of the isotope itself.
  params.addParam<Real>("diffusion_coefficient_base",
                        0.0,
                        "The constant value of the isotope's diffusion "
                        "coefficient.");
  params.addParam<Real>("half_life", 0.0, "The half-life of the isotope (in years).");
  params.addParam<MooseEnum>(
      "half_life_units",
      MooseEnum("seconds minutes hours days years", "years"),
      "The units of the half-lives. Must be consistent for all provided half-lives.");
  params.addParam<std::vector<Real>>("absorption_cross_sections",
                                     "The microscopic absorption cross-sections "
                                     "for the isotope. These must be arranged "
                                     "in increasing order by group.");

  //----------------------------------------------------------------------------
  // Parameters for decay chains.
  params.addParam<std::vector<VariableName>>("decay_parents",
                                             "The variable names of the isotopes "
                                             "which decay into this isotope.");
  params.addParam<std::vector<Real>>("parent_branching_factors",
                                     "The branching factors for the parent "
                                     "isotopes. The order of the branching "
                                     "factors must match the order of "
                                     "decay_parents.");
  params.addParam<std::vector<Real>>("parent_half_lives",
                                     "The half-lives of the parent isotopes "
                                     "(in years). The order of the half-"
                                     "lives must match the order of "
                                     "decay_parents.");

  //----------------------------------------------------------------------------
  // Parameters for neutron activation.
  params.addParam<std::vector<VariableName>>("activation_parents",
                                             "The variable names of the "
                                             "isotopes which are activated to "
                                             "form this isotope.");
  params.addParam<std::vector<Real>>("activation_cross_sections",
                                     "The microscopic activation cross-sections "
                                     "for the isotopes that become this isotope "
                                     "under neutron bombardment. The "
                                     "cross-sections must be ordered by isotope "
                                     "first, neutron flux group second.");

  //----------------------------------------------------------------------------
  // Initial conditions.
  params.addParam<MooseEnum>("ic_type",
                             MooseEnum("constant function file", "constant"),
                             "The type of initial condition to use (if "
                             "multiple are provided). Defaults to constant "
                             "initial conditions.");
  params.addParam<Real>("constant_ic", 0.0, "A constant initial condition for the isotope.");

  //----------------------------------------------------------------------------
  // Parameters required for mass transport and stabilization. These are applied
  // to all isotope kernels.
  params += SetupNuclideSystemAction::validParams();
  params.makeParamNotRequired<MooseEnum>("velocity_type");
  params.makeParamNotRequired<std::vector<VariableName>>("nuclides");

  return params;
}

AddMobileIsotopeAction::AddMobileIsotopeAction(const InputParameters & params)
  : GnatBaseAction(params),
    _hl_units(getParam<MooseEnum>("half_life_units").getEnum<HalfLifeUnits>()),
    _isotope_name(getParam<VariableName>("isotope_name")),
    _decay_parents(getParam<std::vector<VariableName>>("decay_parents")),
    _activation_parents(getParam<std::vector<VariableName>>("activation_parents")),
    _master_isotope_list(getParam<std::vector<VariableName>>("nuclides")),
    _diffusion_coefficient_base(getParam<Real>("diffusion_coefficient_base")),
    _sigma_a(getParam<std::vector<Real>>("absorption_cross_sections")),
    _half_life(0.0),
    _parent_branching_fractions(getParam<std::vector<Real>>("parent_branching_factors")),
    _parent_sigma_act(getParam<std::vector<Real>>("activation_cross_sections")),
    _first_action(true)
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

  // Check to make sure all required isotopes are in the master list.
  {
    auto r = std::find(_master_isotope_list.begin(), _master_isotope_list.end(), _isotope_name);
    if (r == _master_isotope_list.end())
      mooseError("Current isotope " + _isotope_name + " is not in the isotope list.");
  }
  for (const auto & d_parent : _decay_parents)
  {
    auto r = std::find(_master_isotope_list.begin(), _master_isotope_list.end(), d_parent);
    if (r == _master_isotope_list.end())
      mooseError("Parent " + d_parent + " is not in the isotope list.");
  }
  for (const auto & a_parent : _activation_parents)
  {
    auto r = std::find(_master_isotope_list.begin(), _master_isotope_list.end(), a_parent);
    if (r == _master_isotope_list.end())
      mooseError("Parent " + a_parent + " is not in the isotope list.");
  }

  const auto & parent_half_lives = getParam<std::vector<Real>>("parent_half_lives");
  _parent_decay_constants.reserve(parent_half_lives.size());
  switch (_hl_units)
  {
    case HalfLifeUnits::Seconds:
      _half_life = getParam<Real>("half_life");
      for (auto & hl : parent_half_lives)
        _parent_decay_constants.emplace_back(std::log(2.0) / (hl));

      break;

    case HalfLifeUnits::Minutes:
      _half_life = getParam<Real>("half_life") * 60.0;
      for (auto & hl : parent_half_lives)
        _parent_decay_constants.emplace_back(std::log(2.0) / (hl * 60.0));

      break;

    case HalfLifeUnits::Hours:
      _half_life = getParam<Real>("half_life") * 60.0 * 60.0;
      for (auto & hl : parent_half_lives)
        _parent_decay_constants.emplace_back(std::log(2.0) / (hl * 60.0 * 60.0));

      break;

    case HalfLifeUnits::Days:
      _half_life = getParam<Real>("half_life") * 60.0 * 60.0 * 24.0;
      for (auto & hl : parent_half_lives)
        _parent_decay_constants.emplace_back(std::log(2.0) / (hl * 60.0 * 60.0 * 24.0));

      break;

    case HalfLifeUnits::Years:
      _half_life = getParam<Real>("half_life") * 60.0 * 60.0 * 24.0 * 365.0;
      for (auto & hl : parent_half_lives)
        _parent_decay_constants.emplace_back(std::log(2.0) / (hl * 60.0 * 60.0 * 24.0));

      break;

    default:
      break;
  }
}

void
AddMobileIsotopeAction::applyIsotopeParameters(InputParameters & params)
{
  params.set<MooseEnum>("velocity_type") = getParam<MooseEnum>("velocity_type");
  params.set<MooseEnum>("density_type") = getParam<MooseEnum>("density_type");
  params.set<Real>("molar_mass") = getParam<Real>("molar_mass");

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
            params.set<VariableName>("u_var") = getParam<VariableName>("u_var");

            break;

          case ProblemType::Cartesian2D:
            params.set<VariableName>("u_var") = getParam<VariableName>("u_var");
            params.set<VariableName>("v_var") = getParam<VariableName>("v_var");

            break;

          case ProblemType::Cartesian3D:
            params.set<VariableName>("u_var") = getParam<VariableName>("u_var");
            params.set<VariableName>("v_var") = getParam<VariableName>("v_var");
            params.set<VariableName>("w_var") = getParam<VariableName>("w_var");

            break;

          default:
            break;
        }
      }
      else
        params.set<VariableName>("vector_velocity") = getParam<VariableName>("vector_velocity");
      break;
  }
}

void
AddMobileIsotopeAction::addICs()
{
  if (_exec_type != ExecutionType::Transient)
    return;

  switch (static_cast<int>(getParam<MooseEnum>("ic_type")))
  {
    case 0:
      // Add ConstantIC.
      {
        auto params = _factory.getValidParams("ConstantIC");
        params.set<VariableName>("variable") = _isotope_name;

        if (isParamValid("block"))
        {
          params.set<std::vector<SubdomainName>>("block") =
              getParam<std::vector<SubdomainName>>("block");
        }

        params.set<Real>("value") = getParam<Real>("constant_ic");

        _problem->addInitialCondition("ConstantIC", "ConstantIC_" + _isotope_name, params);
        debugOutput("Adding IC ConstantIC for the variable " + _isotope_name + ".");
      } // ConstantIC
      break;

    case 1:
      mooseError("Function ICs are not currently supported.");
      break;

    case 2:
      mooseError("File ICs are not currently supported.");
      break;

    default:
      mooseError("Unknown IC type.");
      break;
  }
}

void
AddMobileIsotopeAction::addKernels()
{
  // Add ADIsotopeActivation.
  if (_activation_parents.size() > 0u)
  {
    auto params = _factory.getValidParams("ADIsotopeActivation");
    params.set<NonlinearVariableName>("variable") = _isotope_name;

    // Apply common isotope parameters.
    applyIsotopeParameters(params);

    // Set the number of groups.
    params.set<unsigned int>("num_groups") = _num_groups;
    // Set the activation cross-sections.
    params.set<std::vector<Real>>("group_activation") = _parent_sigma_act;
    // Apply the parent isotope densities.
    params.set<std::vector<VariableName>>("isotope_densities") = _activation_parents;
    // Apply the scalar fluxes.
    auto & scalar_flux_names = params.set<std::vector<VariableName>>("group_scalar_fluxes");
    scalar_flux_names.reserve(_num_groups);
    for (unsigned int g = 0u; g < _num_groups; ++g)
      scalar_flux_names.emplace_back(_group_flux_moments[g][g * _num_group_moments]);

    if (isParamValid("block"))
    {
      params.set<std::vector<SubdomainName>>("block") =
          getParam<std::vector<SubdomainName>>("block");
    }

    _problem->addKernel("ADIsotopeActivation", "ADIsotopeActivation" + _isotope_name, params);
    debugOutput("    - Adding kernel ADIsotopeActivation for the variable " + _isotope_name + ".");
  } // ADIsotopeActivation

  // Add ADIsotopeAdvection.
  {
    auto params = _factory.getValidParams("ADIsotopeAdvection");
    params.set<NonlinearVariableName>("variable") = _isotope_name;

    // Apply common isotope parameters.
    applyIsotopeParameters(params);

    if (isParamValid("block"))
    {
      params.set<std::vector<SubdomainName>>("block") =
          getParam<std::vector<SubdomainName>>("block");
    }

    _problem->addKernel("ADIsotopeAdvection", "ADIsotopeAdvection" + _isotope_name, params);
    debugOutput("    - Adding kernel ADIsotopeAdvection for the variable " + _isotope_name + ".");
  } // ADIsotopeAdvection

  // Add ADIsotopeDecaySink.
  {
    auto params = _factory.getValidParams("ADIsotopeDecaySink");
    params.set<NonlinearVariableName>("variable") = _isotope_name;

    // Apply common isotope parameters.
    applyIsotopeParameters(params);

    // Decay constant.
    params.set<Real>("decay_const") = std::log(2.0) / _half_life;

    if (isParamValid("block"))
    {
      params.set<std::vector<SubdomainName>>("block") =
          getParam<std::vector<SubdomainName>>("block");
    }

    _problem->addKernel("ADIsotopeDecaySink", "ADIsotopeDecaySink" + _isotope_name, params);
    debugOutput("    - Adding kernel ADIsotopeDecaySink for the variable " + _isotope_name + ".");
  } // ADIsotopeDecaySink

  // Add ADIsotopeDecaySource.
  if (_decay_parents.size() > 0u)
  {
    auto params = _factory.getValidParams("ADIsotopeDecaySource");
    params.set<NonlinearVariableName>("variable") = _isotope_name;

    // Apply common isotope parameters.
    applyIsotopeParameters(params);

    // Apply decay constants and branching factors.
    params.set<std::vector<Real>>("decay_constants") = _parent_decay_constants;
    params.set<std::vector<Real>>("branching_factors") = _parent_branching_fractions;
    // Apply the parent isotope densities.
    params.set<std::vector<VariableName>>("isotope_densities") = _decay_parents;

    if (isParamValid("block"))
    {
      params.set<std::vector<SubdomainName>>("block") =
          getParam<std::vector<SubdomainName>>("block");
    }

    _problem->addKernel("ADIsotopeDecaySource", "ADIsotopeDecaySource" + _isotope_name, params);
    debugOutput("    - Adding kernel ADIsotopeDecaySource for the variable " + _isotope_name + ".");
  } // ADIsotopeDecaySource

  // Add ADIsotopeDepletion.
  if (_sigma_a.size() > 0u)
  {
    auto params = _factory.getValidParams("ADIsotopeDepletion");
    params.set<NonlinearVariableName>("variable") = _isotope_name;

    // Apply common isotope parameters.
    applyIsotopeParameters(params);

    // Set the number of groups.
    params.set<unsigned int>("num_groups") = _num_groups;
    // Set the microscopic absorption cross-sections for this isotope.
    params.set<std::vector<Real>>("group_absorption") = _sigma_a;
    // Apply the scalar fluxes.
    auto & scalar_flux_names = params.set<std::vector<VariableName>>("group_scalar_fluxes");
    scalar_flux_names.reserve(_num_groups);
    for (unsigned int g = 0u; g < _num_groups; ++g)
      scalar_flux_names.emplace_back(_group_flux_moments[g][g * _num_group_moments]);

    if (isParamValid("block"))
    {
      params.set<std::vector<SubdomainName>>("block") =
          getParam<std::vector<SubdomainName>>("block");
    }

    _problem->addKernel("ADIsotopeDepletion", "ADIsotopeDepletion" + _isotope_name, params);
    debugOutput("    - Adding kernel ADIsotopeDepletion for the variable " + _isotope_name + ".");
  } // ADIsotopeDepletion

  // Add ADIsotopeDiffusion.
  {
    auto params = _factory.getValidParams("ADIsotopeDiffusion");
    params.set<NonlinearVariableName>("variable") = _isotope_name;

    // Apply common isotope parameters.
    applyIsotopeParameters(params);

    if (isParamValid("block"))
    {
      params.set<std::vector<SubdomainName>>("block") =
          getParam<std::vector<SubdomainName>>("block");
    }

    _problem->addKernel("ADIsotopeDiffusion", "ADIsotopeDiffusion" + _isotope_name, params);
    debugOutput("    - Adding kernel ADIsotopeDiffusion for the variable " + _isotope_name + ".");
  } // ADIsotopeDiffusion

  // Add ADIsotopeTimeDerivative

  if (_exec_type == ExecutionType::Transient)
  {
    auto params = _factory.getValidParams("ADIsotopeTimeDerivative");
    params.set<NonlinearVariableName>("variable") = _isotope_name;

    applyIsotopeParameters(params);

    if (isParamValid("block"))
    {
      params.set<std::vector<SubdomainName>>("block") =
          getParam<std::vector<SubdomainName>>("block");
    }

    _problem->addKernel(
        "ADIsotopeTimeDerivative", "ADIsotopeTimeDerivative" + _isotope_name, params);
    debugOutput("    - Adding kernel ADIsotopeTimeDerivative for the variable " + _isotope_name +
                ".");
  } // ADIsotopeTimeDerivative
}

void
AddMobileIsotopeAction::addMaterials()
{
  { // Add AutoIsotopeMaterial.
    auto params = _factory.getValidParams("AutoIsotopeMaterial");
    params.set<NonlinearVariableName>("isotope_name") = _isotope_name;

    // TODO: Fluid dependant diffusion coefficient.
    params.set<MooseEnum>("diffusion_coefficient_type") = MooseEnum("constant", "constant");

    params.set<Real>("diffusion_coefficient_base") = _diffusion_coefficient_base;

    _problem->addMaterial("AutoIsotopeMaterial", "AutoIsotopeMaterial" + _isotope_name, params);
    debugOutput("    - Adding material AutoIsotopeMaterial for the variable " + _isotope_name +
                ".");
  } // AutoIsotopeMaterial
}

void
AddMobileIsotopeAction::act()
{
  if (_first_action)
  {
    const std::string init_mass =
        "Initializing a mass transport system for species " + _isotope_name + ": ";
    debugOutput(init_mass, init_mass);

    initializeBase();

    _first_action = false;
  }

  if (_current_task == "add_variable")
  {
    debugOutput("  - Adding variables...");
    addVariable(_isotope_name);
  }

  if (_current_task == "add_ic")
  {
    debugOutput("  - Adding initial conditions...");
    addICs();
  }

  if (_current_task == "add_kernel")
  {
    debugOutput("  - Adding kernels...");
    addKernels();
  }

  if (_current_task == "add_material")
  {
    debugOutput("  - Adding materials...");
    addMaterials();
  }
}
