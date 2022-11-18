#include "NuclideSystemAction.h"

#include <filesystem>
#include <algorithm>

#include "AddVariableAction.h"
#include "SetupNuclideSystemAction.h"
#include "ADIsotopeBase.h"
#include "InputParameterWarehouse.h"

#include "SimpleCSVReader.h"

registerMooseAction("GnatApp", NuclideSystemAction, "add_variable");
registerMooseAction("GnatApp", NuclideSystemAction, "add_ic");
registerMooseAction("GnatApp", NuclideSystemAction, "add_kernel");
registerMooseAction("GnatApp", NuclideSystemAction, "add_material");

InputParameters
NuclideSystemAction::validParams()
{
  auto params = GnatBaseAction::validParams();
  params.addClassDescription(
      "This action adds variables which represents the density of all mobile isotopes "
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

  //----------------------------------------------------------------------------
  // Cross-section file information.
  params.addRequiredParam<std::string>("xs_file_name", "The file to extract cross-sections from.");
  params.addRequiredParam<std::string>("xs_source_material_id",
                                       "The material ID used by the cross-section source to "
                                       "identify cross-sections for different geometric regions. "
                                       "This parameter will be cast to a different type depending "
                                       "on the cross-section source.");
  params.addParam<MooseEnum>("xs_source",
                             MooseEnum("detect gnat openmc", "detect"),
                             "The source which generated the cross-sections. Used for "
                             "parsing different cross-section formats. If set to 'detect' the "
                             "material will attempt to determine the source and parse "
                             "accordingly.");
  params.addParam<MooseEnum>(
      "xs_type",
      MooseEnum("micro macro", "micro"),
      "The type of cross-section in the cross-section file. Can either be microscopic "
      "cross-sections or macroscopic cross-sections. Microscopic is preferred and the default.");

  //----------------------------------------------------------------------------
  // Nuclide system file information.
  params.addRequiredParam<std::string>("nuclide_prop_file_name",
                                       "The file to extract nuclide information from.");
  params.addParam<MooseEnum>(
      "half_life_units",
      MooseEnum("seconds minutes hours days years", "years"),
      "The units of the half-lives. Must be consistent for all provided half-lives.");

  //----------------------------------------------------------------------------
  // Initial conditions and density.
  params.addRequiredParam<Real>("density",
                                "The initial density of the transient nuclide system. The units of "
                                "density must match the units of the mesh.");
  params.addParam<MooseEnum>("ic_type",
                             MooseEnum("constant function file", "constant"),
                             "The type of initial condition to use (if "
                             "multiple are provided). Defaults to constant "
                             "initial conditions. The density multiplied by the weight fraction of "
                             "the nuclide is the constant initial condition.");

  //----------------------------------------------------------------------------
  // Parameters required for mass transport and stabilization. These are applied
  // to all isotope kernels.
  params += SetupNuclideSystemAction::validParams();
  params.makeParamNotRequired<MooseEnum>("velocity_type");
  params.makeParamNotRequired<std::vector<VariableName>>("nuclides");

  return params;
}

NuclideSystemAction::NuclideSystemAction(const InputParameters & parameters)
  : GnatBaseAction(parameters),
    _density(getParam<Real>("density")),
    _xs_source(getParam<MooseEnum>("xs_source").getEnum<CrossSectionSource>()),
    _xs_type(getParam<MooseEnum>("xs_type").getEnum<CrossSectionType>()),
    _xs_file_name(getParam<std::string>("xs_file_name")),
    _xs_source_material_id(getParam<std::string>("xs_source_material_id")),
    _hl_units(getParam<MooseEnum>("half_life_units").getEnum<HalfLifeUnits>()),
    _nuclide_prop_file_name(getParam<std::string>("nuclide_prop_file_name")),
    _first_action(true)
{
  if (_xs_type == CrossSectionType::Macro)
    mooseError("Parsing of macroscopic cross-sections for mobile nuclides is currently not "
               "supported. Please provide microscopic cross-sections.");

  // Check if a container block exists with isotope parameters. If yes, apply them.
  // FIX THIS: The most janky way to fix this breaking MOOSE change.
  auto isotope_system_actions = _awh.getActions<SetupNuclideSystemAction>();
  if (isotope_system_actions.size() == 1)
  {
    const auto & params = _app.getInputParameterWarehouse().getInputParameters();
    InputParameters & pars(*(params.find(uniqueActionName())->second.get()));
    pars.applyParameters(isotope_system_actions[0]->parameters());
  }

  // Setup the nuclide name list.
  for (const auto & name : getParam<std::vector<VariableName>>("nuclides"))
    _nuclide_properties.try_emplace(name);

  // Parse the nuclide system information and required cross-sections.
  parseNuclideSystem();
  parseCrossSections();

  // Convert all cross-sections to effective microscopic cross-sections.
  switch (_xs_type)
  {
    case CrossSectionType::Micro:
      switch (_xs_source)
      {
        case CrossSectionSource::OpenMC:
          // OpenMC stores microscopic cross-sections in barns. This converts them to cm^2
          for (auto & [nuclide, properties] : _nuclide_properties)
          {
            for (unsigned int g = 0u; g < properties._sigma_a_g.size(); ++g)
              properties._sigma_a_g[g] *= (std::pow(10.0, -24.0));
          }
          break;

        default:
          break;
      }
      break;

    default:
      break;
  }

  // Activation parents.
  for (auto & [nuclide, properties] : _nuclide_properties)
  {
    for (const auto & activation_parent : properties._activation_parents)
      properties._sigma_a_g_parent = _nuclide_properties[activation_parent]._sigma_a_g;
  }

  // Assign the activation and decay parent properties to each nuclide.
  // Decay parents.
  switch (_hl_units)
  {
    case HalfLifeUnits::Seconds:
      for (auto & [nuclide, properties] : _nuclide_properties)
      {
        for (const auto & decay_parent : properties._decay_parents)
        {
          if (_nuclide_properties[decay_parent]._half_life < 0.0)
            continue;

          properties._parent_decay_constants.emplace_back(
              std::log(2.0) / (_nuclide_properties[decay_parent]._half_life));
        }
      }
      break;

    case HalfLifeUnits::Minutes:
      for (auto & [nuclide, properties] : _nuclide_properties)
      {
        for (const auto & decay_parent : properties._decay_parents)
        {
          if (_nuclide_properties[decay_parent]._half_life < 0.0)
            continue;

          properties._parent_decay_constants.emplace_back(
              std::log(2.0) / (_nuclide_properties[decay_parent]._half_life * 60.0));
        }
      }
      break;

    case HalfLifeUnits::Hours:
      for (auto & [nuclide, properties] : _nuclide_properties)
      {
        for (const auto & decay_parent : properties._decay_parents)
        {
          if (_nuclide_properties[decay_parent]._half_life < 0.0)
            continue;

          properties._parent_decay_constants.emplace_back(
              std::log(2.0) / (_nuclide_properties[decay_parent]._half_life * 60.0 * 60.0));
        }
      }
      break;

    case HalfLifeUnits::Days:
      for (auto & [nuclide, properties] : _nuclide_properties)
      {
        for (const auto & decay_parent : properties._decay_parents)
        {
          if (_nuclide_properties[decay_parent]._half_life < 0.0)
            continue;

          properties._parent_decay_constants.emplace_back(
              std::log(2.0) / (_nuclide_properties[decay_parent]._half_life * 60.0 * 60.0 * 24.0));
        }
      }
      break;

    case HalfLifeUnits::Years:
      for (auto & [nuclide, properties] : _nuclide_properties)
      {
        for (const auto & decay_parent : properties._decay_parents)
        {
          if (_nuclide_properties[decay_parent]._half_life < 0.0)
            continue;

          properties._parent_decay_constants.emplace_back(
              std::log(2.0) / (_nuclide_properties[decay_parent]._half_life * 60.0 * 60.0 * 24.0));
        }
      }
      break;

    default:
      break;
  }
}

void
NuclideSystemAction::applyIsotopeParameters(InputParameters & params)
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
NuclideSystemAction::addICs(const std::string & var_name, const NuclideProperties & properties)
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

        params.set<Real>("value") = properties._weight_fraction * _density;

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
NuclideSystemAction::addKernels(const std::string & var_name, const NuclideProperties & properties)
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
  if (properties._half_life > 0.0)
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
NuclideSystemAction::addMaterials(const std::string & var_name,
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
NuclideSystemAction::parseNuclideSystem()
{
  // Open the nuclide system file to parse properties.
  std::ifstream descriptor_file(_nuclide_prop_file_name);
  std::string line;
  if (descriptor_file.is_open())
  {
    bool is_block = false;
    std::string nuclide;
    while (std::getline(descriptor_file, line))
    {
      // Skip over lines with a '#' character as they're comments.
      if (line.size() > 0u)
        if (line[0u] == '#')
          continue;

      if (!is_block)
      {
        auto pos = line.find("Begin:");
        if (pos != std::string::npos)
        {
          // We found a block. Parse the nuclide name and set the flag so we acquire the properties.
          line.erase(std::remove(line.begin(), line.end(), ' '), line.end());
          nuclide = line.substr(pos + std::string("Begin:").size());
          is_block = true;
        }
      }
      else if (is_block)
      {
        // We are inside a block. Parse the nuclide properties and check to see if the block ends.
        line.erase(std::remove(line.begin(), line.end(), ' '), line.end());

        // Parse the diffusion coefficient.
        auto pos = line.find("DiffusionCoefficient:");
        if (pos != std::string::npos)
        {
          if (_nuclide_properties.count(nuclide) > 0u)
          {
            _nuclide_properties[nuclide]._diffusion =
                std::stod(line.substr(pos + std::string("DiffusionCoefficient:").size()));
          }
        }

        // Parse the half-life.
        pos = line.find("HalfLife:");
        if (pos != std::string::npos)
        {
          if (_nuclide_properties.count(nuclide) > 0u)
          {
            _nuclide_properties[nuclide]._half_life =
                std::stod(line.substr(pos + std::string("HalfLife:").size()));
          }
        }

        // Parse the weight fraction.
        pos = line.find("WeightFraction:");
        if (pos != std::string::npos)
        {
          if (_nuclide_properties.count(nuclide) > 0u)
          {
            _nuclide_properties[nuclide]._weight_fraction =
                std::stod(line.substr(pos + std::string("WeightFraction:").size()));
          }
        }

        // Parse the decay parent branching factor.
        pos = line.find("DecayParentBranching:");
        if (pos != std::string::npos)
        {
          if (_nuclide_properties.count(nuclide) > 0u)
          {
            auto current_delim_pos = line.find(",");
            std::size_t previous_delim_pos = 0u;
            std::size_t offset = pos + std::string("DecayParentBranching:").size();
            if (current_delim_pos == std::string::npos)
            {
              _nuclide_properties[nuclide]._parent_branching_fractions.emplace_back(
                  std::stod(line.substr(previous_delim_pos + offset)));
            }
            else
            {
              while (current_delim_pos != std::string::npos)
              {
                _nuclide_properties[nuclide]._parent_branching_fractions.emplace_back(
                    std::stod(line.substr(previous_delim_pos + offset,
                                          current_delim_pos - (previous_delim_pos + offset))));
                previous_delim_pos = current_delim_pos;
                current_delim_pos = line.find(",", previous_delim_pos + 1);
                offset = 1u;
              }
            }
          }
        }

        // Parse the decay parents.
        pos = line.find("DecayParents:");
        if (pos != std::string::npos)
        {
          if (_nuclide_properties.count(nuclide) > 0u)
          {
            auto current_delim_pos = line.find(",");
            std::size_t previous_delim_pos = 0u;
            std::size_t offset = pos + std::string("DecayParents:").size();
            if (current_delim_pos == std::string::npos)
            {
              _nuclide_properties[nuclide]._decay_parents.emplace_back(
                  line.substr(previous_delim_pos + offset));
            }
            else
            {
              while (current_delim_pos != std::string::npos)
              {
                _nuclide_properties[nuclide]._decay_parents.emplace_back(
                    line.substr(previous_delim_pos + offset,
                                current_delim_pos - (previous_delim_pos + offset)));
                previous_delim_pos = current_delim_pos;
                current_delim_pos = line.find(",", previous_delim_pos + 1);
                offset = 1u;
              }
            }
          }
        }

        // Parse the activation parents.
        pos = line.find("ActivationParents:");
        if (pos != std::string::npos)
        {
          if (_nuclide_properties.count(nuclide) > 0u)
          {
            auto current_delim_pos = line.find(",");
            std::size_t previous_delim_pos = 0u;
            std::size_t offset = pos + std::string("ActivationParents:").size();
            if (current_delim_pos == std::string::npos)
            {
              _nuclide_properties[nuclide]._activation_parents.emplace_back(
                  line.substr(previous_delim_pos + offset));
            }
            else
            {
              while (current_delim_pos != std::string::npos)
              {
                _nuclide_properties[nuclide]._activation_parents.emplace_back(
                    line.substr(previous_delim_pos + offset,
                                current_delim_pos - (previous_delim_pos + offset)));
                previous_delim_pos = current_delim_pos;
                current_delim_pos = line.find(",", previous_delim_pos + 1);
                offset = 1u;
              }
            }
          }
        }

        // Case where the block ends.
        pos = line.find(std::string("End:") + nuclide);
        if (pos != std::string::npos)
        {
          // We found the end of a block.
          nuclide = "";
          is_block = false;
        }
      }
    }

    descriptor_file.close();
  }
  else
    mooseError("Failed to open file " + _nuclide_prop_file_name + " to read material properties.");
}

void
NuclideSystemAction::parseCrossSections()
{
  // Open the cross-section descriptor file to figure out where the cross-sections are stored and
  // what the files are named.
  std::ifstream descriptor_file(_xs_file_name);
  std::string line;
  if (descriptor_file.is_open())
  {
    std::getline(descriptor_file, line);
    if (_xs_source == CrossSectionSource::Detect)
    {
      // First line should always be the cross-section source. Use that to determine the
      // properties.
      if (line == "openmc")
        _xs_source = CrossSectionSource::OpenMC;
      else if (line == "gnat")
        _xs_source = CrossSectionSource::Gnat;
      else
        mooseError("Unknown cross-section source.");
    }

    // Read each line. Assume that anything after the keyword is the relative path of the
    // specific cross-section file to the descriptor file.
    while (std::getline(descriptor_file, line))
    {
      auto pos = line.find("SigmaA: ");
      if (pos != std::string::npos)
      {
        parseProperty(PropertyType::SigmaA, line.substr(pos + std::string("SigmaA: ").size()));
        continue;
      }
    }

    descriptor_file.close();
  }
  else
    mooseError("Failed to open file " + _xs_file_name + " to read material properties.");
}

void
NuclideSystemAction::parseProperty(const PropertyType & type, const std::string & property_file)
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
NuclideSystemAction::parseGnatProperty(const PropertyType & type, const std::string & property_file)
{
  mooseError("Parsing of cross-sections generated by Gnat is currently not supported. " +
             property_file);

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
NuclideSystemAction::parseOpenMCProperty(const PropertyType & /*type*/,
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
    {
      std::string first_nuclide = reader.getEntry("nuclide", static_cast<unsigned int>(min_row));
      for (unsigned int row = static_cast<unsigned int>(min_row);
           row <= static_cast<unsigned int>(max_row);
           ++row)
      {
        num_groups = std::max(num_groups, std::stoi(reader.getEntry("group in", row)));

        if (first_nuclide == reader.getEntry("nuclide", row) &&
            row != static_cast<unsigned int>(min_row) && first_nuclide != "")
          first_nuclide = "";
      }
    }

    // Loop over all of the rows of data and parse accordingly.
    // OpenMC orders energy groups from the largest to smallest. As an example: group 1 would be
    // fast, group 2 would be thermal.
    {
      auto & values = reader.getColumn("mean");
      for (unsigned int row = static_cast<unsigned int>(min_row);
           row <= static_cast<unsigned int>(max_row);
           ++row)
      {
        if (_nuclide_properties.count(reader.getEntry("nuclide", row)) > 0u)
        {
          _nuclide_properties[reader.getEntry("nuclide", row)]._sigma_a_g.emplace_back(
              std::stod(values[row]));
        }
      }
    }
  }
  else
    mooseError("Failed to open cross-section file: " + dir.string());
}

void
NuclideSystemAction::act()
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
      addICs(nuclide, properties);
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
