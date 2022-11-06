#include "AddFileIsotopesAction.h"

#include <filesystem>

#include "AddVariableAction.h"
#include "SetupIsotopeSystemAction.h"
#include "ADIsotopeBase.h"
#include "InputParameterWarehouse.h"

#include "SimpleCSVReader.h"

registerMooseAction("GnatApp", AddFileIsotopesAction, "add_variable");
registerMooseAction("GnatApp", AddFileIsotopesAction, "add_ic");
registerMooseAction("GnatApp", AddFileIsotopesAction, "add_kernel");
registerMooseAction("GnatApp", AddFileIsotopesAction, "add_material");

InputParameters
AddFileIsotopesAction::validParams()
{
  auto params = GnatBaseAction::validParams();
  params.addClassDescription(
      "This action adds variables which represents the mass/number density of all mobile isotopes "
      "contained in the provided cross-section files, alongside the kernels required to describe "
      "the time evolution of those isotopes.");

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
  params.addRequiredParam<std::vector<VariableName>>(
      "isotope_names", "The identifying names of the isotopes to pull from the provided files.");

  //----------------------------------------------------------------------------
  // Property file information.
  params.addRequiredParam<std::string>("xs_file_name", "The file to extract cross-sections from.");
  params.addRequiredParam<std::string>("xs_source_material_id",
                                       "The material ID used by the cross-section source to "
                                       "identify cross-sections for different geometric regions. "
                                       "This parameter will be cast to a different type depending "
                                       "on the cross-section source.");
  params.addParam<MooseEnum>("cross_section_source",
                             MooseEnum("detect gnat openmc", "detect"),
                             "The source which generated the cross-sections. Used for "
                             "parsing different cross-section formats. If set to 'detect' the "
                             "material will attempt to determine the source and parse "
                             "accordingly.");
  params.addRequiredParam<std::string>("nuclide_prop_file_name",
                                       "The file to extract nuclide information from.");
  params.addParam<MooseEnum>(
      "half_life_units",
      MooseEnum("seconds minutes hours days years", "years"),
      "The units of the half-lives. Must be consistent for all provided half-lives.");

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
  params += SetupIsotopeSystemAction::validParams();
  params.makeParamNotRequired<MooseEnum>("velocity_type");
  params.makeParamNotRequired<std::vector<VariableName>>("isotopes");

  return params;
}

AddFileIsotopesAction::AddFileIsotopesAction(const InputParameters & parameters)
  : GnatBaseAction(parameters),
    _xs_file_name(getParam<std::string>("xs_file_name")),
    _xs_source_material_id(getParam<std::string>("xs_source_material_id")),
    _xs_source(getParam<MooseEnum>("cross_section_source").getEnum<CrossSectionSource>()),
    _nuclide_prop_file_name(getParam<std::string>("nuclide_prop_file_name")),
    _hl_units(getParam<MooseEnum>("half_life_units").getEnum<HalfLifeUnits>()),
    _first_action(true)
{
  // Check if a container block exists with isotope parameters. If yes, apply them.
  // FIX THIS: The most janky way to fix this breaking MOOSE change.
  auto isotope_system_actions = _awh.getActions<SetupIsotopeSystemAction>();
  if (isotope_system_actions.size() == 1)
  {
    const auto & params = _app.getInputParameterWarehouse().getInputParameters();
    InputParameters & pars(*(params.find(uniqueActionName())->second.get()));
    pars.applyParameters(isotope_system_actions[0]->parameters());
  }

  // Parse cross-sections first.

  // Nuclide information second.
}

void
AddFileIsotopesAction::applyIsotopeParameters(InputParameters & params)
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
AddFileIsotopesAction::addICs(const std::string & var_name)
{
  if (_exec_type != ExecutionType::Transient)
    return;

  switch (static_cast<int>(getParam<MooseEnum>("ic_type")))
  {
    case 0:
      // Add ConstantIC.
      {
        auto params = _factory.getValidParams("ConstantIC");
        params.set<VariableName>("variable") = var_name;

        if (isParamValid("block"))
        {
          params.set<std::vector<SubdomainName>>("block") =
              getParam<std::vector<SubdomainName>>("block");
        }

        params.set<Real>("value") = getParam<Real>("constant_ic");

        _problem->addInitialCondition("ConstantIC", "ConstantIC_" + var_name, params);
        debugOutput("Adding IC ConstantIC for the variable " + var_name + ".");
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
AddFileIsotopesAction::addKernels(const std::string & var_name,
                                  const NuclideProperties & properties)
{
  // Add ADIsotopeActivation.
  if (properties._activation_parents.size() > 0u)
  {
    auto params = _factory.getValidParams("ADIsotopeActivation");
    params.set<NonlinearVariableName>("variable") = var_name;

    // Apply common isotope parameters.
    applyIsotopeParameters(params);

    // Set the number of groups.
    params.set<unsigned int>("num_groups") = _num_groups;
    // Set the activation cross-sections.
    params.set<std::vector<Real>>("group_activation") = properties._sigma_a_g_parent;
    // Apply the parent isotope densities.
    params.set<std::vector<VariableName>>("isotope_densities") = properties._activation_parents;
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

    _problem->addKernel("ADIsotopeActivation", "ADIsotopeActivation" + var_name, params);
    debugOutput("    - Adding kernel ADIsotopeActivation for the variable " + var_name + ".");
  } // ADIsotopeActivation

  // Add ADIsotopeAdvection.
  {
    auto params = _factory.getValidParams("ADIsotopeAdvection");
    params.set<NonlinearVariableName>("variable") = var_name;

    // Apply common isotope parameters.
    applyIsotopeParameters(params);

    if (isParamValid("block"))
    {
      params.set<std::vector<SubdomainName>>("block") =
          getParam<std::vector<SubdomainName>>("block");
    }

    _problem->addKernel("ADIsotopeAdvection", "ADIsotopeAdvection" + var_name, params);
    debugOutput("    - Adding kernel ADIsotopeAdvection for the variable " + var_name + ".");
  } // ADIsotopeAdvection

  // Add ADIsotopeDecaySink.
  {
    auto params = _factory.getValidParams("ADIsotopeDecaySink");
    params.set<NonlinearVariableName>("variable") = var_name;

    // Apply common isotope parameters.
    applyIsotopeParameters(params);

    // Decay constant.
    params.set<Real>("decay_const") = std::log(2.0) / properties._half_life;

    if (isParamValid("block"))
    {
      params.set<std::vector<SubdomainName>>("block") =
          getParam<std::vector<SubdomainName>>("block");
    }

    _problem->addKernel("ADIsotopeDecaySink", "ADIsotopeDecaySink" + var_name, params);
    debugOutput("    - Adding kernel ADIsotopeDecaySink for the variable " + var_name + ".");
  } // ADIsotopeDecaySink

  // Add ADIsotopeDecaySource.
  if (properties._decay_parents.size() > 0u)
  {
    auto params = _factory.getValidParams("ADIsotopeDecaySource");
    params.set<NonlinearVariableName>("variable") = var_name;

    // Apply common isotope parameters.
    applyIsotopeParameters(params);

    // Apply decay constants and branching factors.
    params.set<std::vector<Real>>("decay_constants") = properties._parent_decay_constants;
    params.set<std::vector<Real>>("branching_factors") = properties._parent_branching_fractions;
    // Apply the parent isotope densities.
    params.set<std::vector<VariableName>>("isotope_densities") = properties._decay_parents;

    if (isParamValid("block"))
    {
      params.set<std::vector<SubdomainName>>("block") =
          getParam<std::vector<SubdomainName>>("block");
    }

    _problem->addKernel("ADIsotopeDecaySource", "ADIsotopeDecaySource" + var_name, params);
    debugOutput("    - Adding kernel ADIsotopeDecaySource for the variable " + var_name + ".");
  } // ADIsotopeDecaySource

  // Add ADIsotopeDepletion.
  if (properties._sigma_a_g.size() > 0u)
  {
    auto params = _factory.getValidParams("ADIsotopeDepletion");
    params.set<NonlinearVariableName>("variable") = var_name;

    // Apply common isotope parameters.
    applyIsotopeParameters(params);

    // Set the number of groups.
    params.set<unsigned int>("num_groups") = _num_groups;
    // Set the microscopic absorption cross-sections for this isotope.
    params.set<std::vector<Real>>("group_absorption") = properties._sigma_a_g;
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

    _problem->addKernel("ADIsotopeDepletion", "ADIsotopeDepletion" + var_name, params);
    debugOutput("    - Adding kernel ADIsotopeDepletion for the variable " + var_name + ".");
  } // ADIsotopeDepletion

  // Add ADIsotopeDiffusion.
  {
    auto params = _factory.getValidParams("ADIsotopeDiffusion");
    params.set<NonlinearVariableName>("variable") = var_name;

    // Apply common isotope parameters.
    applyIsotopeParameters(params);

    if (isParamValid("block"))
    {
      params.set<std::vector<SubdomainName>>("block") =
          getParam<std::vector<SubdomainName>>("block");
    }

    _problem->addKernel("ADIsotopeDiffusion", "ADIsotopeDiffusion" + var_name, params);
    debugOutput("    - Adding kernel ADIsotopeDiffusion for the variable " + var_name + ".");
  } // ADIsotopeDiffusion

  // Add ADIsotopeTimeDerivative

  if (_exec_type == ExecutionType::Transient)
  {
    auto params = _factory.getValidParams("ADIsotopeTimeDerivative");
    params.set<NonlinearVariableName>("variable") = var_name;

    applyIsotopeParameters(params);

    if (isParamValid("block"))
    {
      params.set<std::vector<SubdomainName>>("block") =
          getParam<std::vector<SubdomainName>>("block");
    }

    _problem->addKernel("ADIsotopeTimeDerivative", "ADIsotopeTimeDerivative" + var_name, params);
    debugOutput("    - Adding kernel ADIsotopeTimeDerivative for the variable " + var_name + ".");
  } // ADIsotopeTimeDerivative
}

void
AddFileIsotopesAction::addMaterials(const std::string & var_name,
                                    const NuclideProperties & properties)
{
  { // Add AutoIsotopeMaterial.
    auto params = _factory.getValidParams("AutoIsotopeMaterial");
    params.set<NonlinearVariableName>("isotope_name") = var_name;

    // TODO: Fluid dependant diffusion coefficient.
    params.set<MooseEnum>("diffusion_coefficient_type") = MooseEnum("constant", "constant");

    params.set<Real>("diffusion_coefficient_base") = properties._diffusion;

    _problem->addMaterial("AutoIsotopeMaterial", "AutoIsotopeMaterial" + var_name, params);
    debugOutput("    - Adding material AutoIsotopeMaterial for the variable " + var_name + ".");
  } // AutoIsotopeMaterial
}

void
AddFileIsotopesAction::parseProperty(const PropertyType & type, const std::string & property_file)
{
  switch (_xs_source)
  {
    case CrossSectionSource::Gnat:
      parseGnatProperty(type, property_file);
      break;

    case CrossSectionSource::OpenMC:
      parseOpenMCProperty(type, property_file);
      break;

    default:
      break;
  }
}

void
AddFileIsotopesAction::parseGnatProperty(const PropertyType & type,
                                         const std::string & property_file)
{
  mooseError("Parsing of cross-sections generated by Gnat is currently not supported.");

  switch (type)
  {
    case PropertyType::InvVelocity:
      break;

    case PropertyType::SigmaR:
      break;

    case PropertyType::SigmaA:
      break;

    case PropertyType::SigmaS:
      break;

    case PropertyType::SigmaSMatrix:
      break;

    default:
      break;
  }
}

void
AddFileIsotopesAction::parseOpenMCProperty(const PropertyType & type,
                                           const std::string & property_file)
{
  // Compute the actual path of the property file.
  auto dir = std::filesystem::path(_xs_file_name).parent_path();
  dir /= std::filesystem::path(property_file);

  // Parse!
  SimpleCSVReader reader(dir.string());
  if (reader.read())
  {
    // Loop over the cell ids to find ones that match the provided id in the material.
    int min_row = -1;
    int max_row = -1;
    {
      auto & cell_ids = reader.getColumn("cell");
      // Start by finding the minimum row.
      for (unsigned int row = 0u; row < cell_ids.size(); ++row)
      {
        if (cell_ids[row] == _xs_source_material_id)
        {
          min_row = static_cast<int>(row);
          break;
        }
      }

      if (min_row == -1)
        mooseError("OpenMC cell with label " + _xs_source_material_id + " does not exist in " +
                   property_file);

      // Next find the maximum row (taking advantage of how OpenMC groups cells together).
      for (unsigned int row = static_cast<unsigned int>(min_row) + 1; row < cell_ids.size(); ++row)
      {
        if (cell_ids[row - 1u] == _xs_source_material_id && cell_ids[row] != _xs_source_material_id)
        {
          max_row = static_cast<int>(row - 1u);
          break;
        }
      }

      if (min_row != -1 && max_row == -1)
        max_row = static_cast<int>(cell_ids.size() - 1u);
    }

    // Loop over the data to determine the number of isotopes.
    int num_groups = -1;
    unsigned int num_isotopes = 0;
    {
      std::string first_nuclide = reader.getEntry("nuclide", static_cast<unsigned int>(min_row));
      for (unsigned int row = static_cast<unsigned int>(min_row);
           row <= static_cast<unsigned int>(max_row);
           ++row)
      {
        num_groups = std::max(num_groups, std::stoi(reader.getEntry("group in", row)));

        if (first_nuclide == reader.getEntry("nuclide", row) &&
            row != static_cast<unsigned int>(min_row) && first_nuclide != "")
        {
          num_isotopes = row - static_cast<unsigned int>(min_row);
          first_nuclide = "";
        }
      }
    }

    // Loop over the data and add nuclides to the unordered map of data.
    {
      auto & nuclides = reader.getColumn("nuclide");
      for (unsigned int row = static_cast<unsigned int>(min_row);
           row < static_cast<unsigned int>(min_row) + num_isotopes;
           ++row)
        _nuclide_properties.try_emplace(nuclides[row]);
    }

    // Loop over all of the rows of data and parse accordingly.
    // OpenMC orders energy groups from the largest to smallest. As an example: group 1 would be
    // fast, group 2 would be thermal.
    {
      auto & values = reader.getColumn("mean");
      for (unsigned int row = static_cast<unsigned int>(min_row);
           row <= static_cast<unsigned int>(max_row);
           ++row)
        _nuclide_properties[reader.getEntry("nuclide", row)]._sigma_a_g.emplace_back(
            std::stod(values[row]));
    }
  }
  else
    mooseError("Failed to open cross-section file: " + dir.string());
}

void
AddFileIsotopesAction::act()
{
  if (_first_action)
  {
    debugOutput("Initializing a mass transport system: ", "Initializing a mass transport system: ");

    initializeBase();

    _first_action = false;
  }

  for (const auto & [nuclide, properties] : _nuclide_properties)
  {
    if (_current_task == "add_variable")
    {
      debugOutput("  - Adding variables...");
      addVariable(nuclide);
    }

    if (_current_task == "add_ic")
    {
      debugOutput("  - Adding initial conditions...");
      addICs(nuclide);
    }

    if (_current_task == "add_kernel")
    {
      debugOutput("  - Adding kernels...");
      addKernels(nuclide, properties);
    }

    if (_current_task == "add_material")
    {
      debugOutput("  - Adding materials...");
      addMaterials(nuclide, properties);
    }
  }
}
