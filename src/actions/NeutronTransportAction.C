#include "NeutronTransportAction.h"

#include "Factory.h"
#include "Parser.h"
#include "NonlinearSystemBase.h"
#include "FEProblemBase.h"
#include "ConsoleStream.h"

#include "Conversion.h"
#include "MooseTypes.h"
#include "FEProblem.h"

#include "AddVariableAction.h"
#include "ADDFEMUpwinding.h"

#include "libmesh/string_to_enum.h"
#include "libmesh/fe_type.h"

registerMooseAction("GnatApp", NeutronTransportAction, "add_variable");
registerMooseAction("GnatApp", NeutronTransportAction, "add_kernel");
registerMooseAction("GnatApp", NeutronTransportAction, "add_dg_kernel");
registerMooseAction("GnatApp", NeutronTransportAction, "add_dirac_kernel");
registerMooseAction("GnatApp", NeutronTransportAction, "add_bc");
registerMooseAction("GnatApp", NeutronTransportAction, "add_ic");
registerMooseAction("GnatApp", NeutronTransportAction, "add_aux_variable");
registerMooseAction("GnatApp", NeutronTransportAction, "add_aux_kernel");
registerMooseAction("GnatApp", NeutronTransportAction, "add_output");

InputParameters
NeutronTransportAction::validParams()
{
  auto params = Action::validParams();
  params.addClassDescription("This action adds all of the required variables, "
                             "kernels, boundary conditions, initial conditions "
                             "and auxiliary systems required to solve "
                             "source-driven multi-group neutron transport "
                             "problems with Gnat.");

  //----------------------------------------------------------------------------
  // Parameters for variables.
  params.addRequiredParam<MooseEnum>("family", AddVariableAction::getNonlinearVariableFamilies(),
                                     "Specifies the family of FE shape functions to "
                                     "use for this variable.");
  params.addRequiredParam<MooseEnum>("order", AddVariableAction::getNonlinearVariableOrders(),
                                     "Specifies the order of the FE shape "
                                     "function to use for this variable "
                                     "(additional orders not listed are "
                                     "allowed).");
  params.addParam<Real>("scaling", 1.0,
                        "Specifies a scaling factor to apply to "
                        "this variable.");
  params.addParam<std::string>("angular_flux_names",
                               "angular_flux",
                               "Variable names for the angular flux. The output "
                               "format for the group angular fluxes will be of "
                               "the form {angular_flux_names}_g_n.");
  params.addParam<std::string>("flux_moment_names",
                               "flux_moment",
                               "Variable names for the moments of the angular "
                               "flux. The output format for the group flux "
                               "moments will be of the form "
                               "{flux_moment_names}_g_l_m.");
  params.addParam<bool>("output_angular_fluxes", false,
                        "Whether the angular flux ordinates should be written "
                        "to the exodus file or not");
  params.addParamNamesToGroup("scaling angular_flux_names "
                              "flux_moment_names output_angular_fluxes",
                              "Variable");

  //----------------------------------------------------------------------------
  // Basic parameters for the neutron transport simulation.
  params.addRequiredRangeCheckedParam<unsigned int>("num_groups",
                                                    "num_groups > 0",
                                                    "The number of spectral "
                                                    "energy groups in the "
                                                    "problem.");
  params.addRequiredParam<MooseEnum>("scheme",
                                     MooseEnum("saaf_cfem upwinding_dfem"),
                                     "The discretization and stabilization "
                                     "scheme for the transport equation.");
  params.addRequiredParam<MooseEnum>("execution_type",
                                     MooseEnum("steady transient"),
                                     "The method of execution for the problem. "
                                     "Options are steady-state source driven "
                                     "problems and transient source problems.");
  params.addParam<unsigned int>("max_anisotropy", 0,
                                "The maximum degree of anisotropy to evaluate. "
                                "Defaults to 0 for isotropic scattering.");
  params.addParam<std::vector<SubdomainName>>("block",
                                              "The list of blocks (ids or "
                                              "names) that this variable will "
                                              "be applied.");
  params.addParamNamesToGroup("max_anisotropy block", "Simulation");

  //----------------------------------------------------------------------------
  // Quadrature parameters.
  params.addRangeCheckedParam<unsigned int>("n_polar", 3, "n_polar > 0",
                                            "Number of Legendre polar "
                                            "quadrature points in a single "
                                            "octant of the unit sphere. "
                                            "Defaults to 3.");
  params.addRangeCheckedParam<unsigned int>("n_azimuthal", 3, "n_azimuthal > 0",
                                            "Number of Chebyshev azimuthal "
                                            "quadrature points in a single "
                                            "octant of the unit sphere. "
                                            "Defaults to 3.");
  params.addParam<MooseEnum>("major_axis", MooseEnum("x y z", "x"),
                             "Major axis of the angular quadrature. Allows the "
                             "polar angular quadrature to align with a cartesian "
                             "axis with minimal heterogeneity. Default is the "
                             "x-axis. This parameter is only applied in 3D "
                             "cartesian problems.");
  params.addParamNamesToGroup("n_polar n_azimuthal major_axis",
                              "Quadrature");

  //----------------------------------------------------------------------------
  // Boundary conditions.
  params.addParam<std::vector<BoundaryName>>("vacuum_boundaries",
                                             "The boundaries to apply vacuum "
                                             "boundary conditions.");
  params.addParam<std::vector<BoundaryName>>("source_boundaries",
                                             "The boundaries to apply incoming "
                                             "flux boundary conditions.");
  params.addParam<std::vector<BoundaryName>>("reflective_boundaries",
                                             "The boundaries to apply reflective "
                                             "boundary conditions.");
  params.addParamNamesToGroup("vacuum_boundaries source_boundaries "
                              "reflective_boundaries", "Boundary Condition");

  //----------------------------------------------------------------------------
  // Isotropic point sources.
  params.addParam<std::vector<Point>>("point_source_locations",
                                      "The locations of all isotropic "
                                      "point sources in the problem "
                                      "space.");
  params.addParam<std::vector<Real>>("point_source_intensities",
                                     "The intensities of all "
                                     "isotropic point sources.");
  params.addParam<std::vector<unsigned int>>("point_source_groups",
                                             "The spectral energy groups each "
                                             "isotropic point source emits "
                                             "into.");
  params.addParamNamesToGroup("point_source_locations point_source_intensities "
                              "point_source_groups", "Isotropic Point Source");

  //----------------------------------------------------------------------------
  // Initial conditions.
  params.addParam<MooseEnum>("ic_type",
                             MooseEnum("constant function file", "constant"),
                             "The type of initial condition to use (if "
                             "multiple are provided). Defaults to constant "
                             "initial conditions.");
  params.addParam<std::vector<Real>>("constant_ic", std::vector<Real>(0.0),
                                     "A constant initial condition for the "
                                     "angular fluxes.");
  params.addParamNamesToGroup("ic_type constant_ic", "Initial Condition");
  // TODO: Init from function.
  // TODO: Init from file.

  //----------------------------------------------------------------------------
  // Parameters for debugging.
  params.addParam<MooseEnum>("debug_verbosity",
                             MooseEnum("level0 level1", "level1"),
                             "How verbose the debug output of the transport "
                             "system should be. level0 is fully verbose. "
                             "level1 outputs less debugging information.");
  params.addParam<bool>("debug_disable_scattering", false,
                        "Debug option to disable scattering evaluation.");
  params.addParam<bool>("debug_disable_source_iteration", true,
                        "Debug option to disable source iteration.");
  params.addParamNamesToGroup("debug_verbosity debug_disable_scattering "
                              "debug_disable_source_iteration", "Debugging");

  params.addParamNamesToGroup("family order num_groups execution_type scheme",
                              "Required");

  return params;
}

NeutronTransportAction::NeutronTransportAction(const InputParameters & params)
  : Action(params)
  , _num_groups(getParam<unsigned int>("num_groups"))
  , _transport_scheme(getParam<MooseEnum>("scheme").getEnum<Scheme>())
  , _exec_type(getParam<MooseEnum>("execution_type").getEnum<ExecutionType>())
  , _n_l(getParam<unsigned int>("n_polar"))
  , _n_c(getParam<unsigned int>("n_azimuthal"))
  , _max_eval_anisotropy(getParam<unsigned int>("max_anisotropy"))
  , _angular_flux_name(getParam<std::string>("angular_flux_names"))
  , _flux_moment_name(getParam<std::string>("flux_moment_names"))
  , _vacuum_side_sets(getParam<std::vector<BoundaryName>>("vacuum_boundaries"))
  , _source_side_sets(getParam<std::vector<BoundaryName>>("source_boundaries"))
  , _reflective_side_sets(getParam<std::vector<BoundaryName>>("reflective_boundaries"))
  , _debug_level(getParam<MooseEnum>("debug_verbosity").getEnum<DebugVerbosity>())
  , _var_init(false)
{ }

void
NeutronTransportAction::addRelationshipManagers(Moose::RelationshipManagerType when_type)
{
  if (_transport_scheme == Scheme::UpwindingDFEM)
  {
    auto params = ADDFEMUpwinding::validParams();
    addRelationshipManagers(when_type, params);
  }
}

void
NeutronTransportAction::act()
{
  // Initialize common members.
  initializeCommon();

  // Act functions for different schemes.
  switch (_transport_scheme)
  {
    case Scheme::SAAFCFEM: actSAAFCFEM(); break;
    case Scheme::UpwindingDFEM: actUpwindDFEM(); break;
    default: break;
  }
}

//------------------------------------------------------------------------------
// Function for debug output.
//------------------------------------------------------------------------------
void
NeutronTransportAction::debugOutput(const std::string & level0,
                                    const std::string & level1)
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
//------------------------------------------------------------------------------

//------------------------------------------------------------------------------
// Function to initialize SN quadrature parameters.
//------------------------------------------------------------------------------
void
NeutronTransportAction::applyQuadratureParameters(InputParameters & params)
{
  // Assign dimensionality.
  MooseEnum dimensionality("1D_cartesian 2D_cartesian 3D_cartesian");
  dimensionality.assign(static_cast<int>(_p_type));
  params.set<MooseEnum>("dimensionality") = dimensionality;

  // Assign major axis.
  params.set<MooseEnum>("major_axis") = getParam<MooseEnum>("major_axis");

  // Assign quadrature rules.
  params.set<unsigned int>("n_l") = 2u * _n_l;
  params.set<unsigned int>("n_c") = 2u * _n_c;
}
//------------------------------------------------------------------------------

//------------------------------------------------------------------------------
// Function to initialize common parameters for all schemes.
//------------------------------------------------------------------------------
void
NeutronTransportAction::initializeCommon()
{
  if (!_var_init)
  {
    // Warn the user about issues with certain schemes.
    switch (_transport_scheme)
    {
      case Scheme::SAAFCFEM:
        mooseWarning("The CGFEM-SAAF scheme currently has issues with negative "
                     "fluxes.");

      break;

      case Scheme::UpwindingDFEM:
        //
        mooseWarning("The DGFEM-Upwinding scheme currently has issues with "
                     "negative fluxes.");

        break;

      default: break;
    }

    debugOutput("Transport System Initialization: ",
                "Transport System Initialization: ");

    switch (_transport_scheme)
    {
      case Scheme::SAAFCFEM:
        debugOutput("  - Scheme: SAAF-CGFEM", "  - Scheme: SAAF-CGFEM");
        break;
      case Scheme::UpwindingDFEM:
        debugOutput("  - Scheme: Upwinding-DGFEM", "  - Scheme: Upwinding-DGFEM");
        break;
    }

    debugOutput("  - Building the angular quadrature set...",
                "  - Building the angular quadrature set...");

    // Setup for all actions the NeutronTransport action performs.
    // Grab all the subdomain IDs that the neutron transport action should be
    // applied to.
    std::vector<SubdomainName> block_names(getParam<std::vector<SubdomainName>>("block"));
    for (const auto & block_name : block_names)
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
        _num_flux_ordinates = 2u * _n_l;
        break;

      case 2u:
        _p_type = ProblemType::Cartesian2D;
        _num_group_moments = (_max_eval_anisotropy + 1u) * (_max_eval_anisotropy + 2u) / 2u;
        _num_flux_ordinates = 4u * _n_l * _n_c;
        break;

      case 3u:
        _p_type = ProblemType::Cartesian3D;
        _num_group_moments = (_max_eval_anisotropy + 1u) * (_max_eval_anisotropy + 1u);
        _num_flux_ordinates = 8u * _n_l * _n_c;
        break;

      default:
        mooseError("Unknown mesh dimensionality.");
        break;
    }

    std::string set_info(std::string("    - Angular quadrature set information:\n") +
                         std::string("      - Polar points per quadrant: ") +
                         Moose::stringify(_n_l) +
                         std::string("\n      - Azimuthal points per quadrant: ") +
                         Moose::stringify(_n_c) +
                         std::string("\n      - Total number of flux ordinates: ") +
                         Moose::stringify(_num_flux_ordinates));
    debugOutput(set_info);

    debugOutput("  - Determining variable names...",
                "  - Determining variable names...");
    for (unsigned int g = 0; g < _num_groups; ++g)
    {
      // Set up variable names for the group angular fluxes.
      _group_angular_fluxes.emplace(g, std::vector<VariableName>());
      _group_angular_fluxes[g].reserve(_num_flux_ordinates);
      for (unsigned int n = 1; n <= _num_flux_ordinates; ++n)
      {
        _group_angular_fluxes[g].emplace_back(_angular_flux_name + "_"
                                             + Moose::stringify(g + 1u)  + "_"
                                             + Moose::stringify(n));
      }

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

    _var_init = true;
    debugOutput("  - Initializing flux groups...", "  - Initializing flux groups...");
  }

  // Loop over all groups.
  for (unsigned int g = 0; g < _num_groups; ++g)
  {
    if (!_var_init)
      debugOutput("  - Initializing flux ordinates...");

    // Loop over all the flux ordinates.
    for (unsigned int n = 0; n < _num_flux_ordinates; ++n)
    {
      const auto & var_name = _group_angular_fluxes[g][n];

      // Add a non-linear variable.
      if (_current_task == "add_variable")
      {
        if (g == 0u && n == 0u)
          debugOutput("    - Adding variables...");

        addVariable(var_name);
      }

      // Add boundary conditions.
      if (_current_task == "add_bc")
      {
        if (g == 0u && n == 0u)
          debugOutput("    - Adding BCs...");

        addBCs(var_name, g, n);
      }

      // Add initial conditions.
      if (_current_task == "add_ic")
      {
        if (g == 0u && n == 0u)
          debugOutput("    - Adding ICs...");

        addICs(var_name, g, n);
      }
    }

    // Loop over all moments and set up auxvariables and auxkernels.
    unsigned int moment_index = 0u;
    for (unsigned int l = 0; l <= _max_eval_anisotropy; ++l)
    {
      for (int m = -1 * static_cast<int>(l); m <= static_cast<int>(l); ++m)
      {
        const auto & var_name = _group_flux_moments[g][moment_index];

        // Add auxvariables.
        if (_current_task == "add_aux_variable")
        {
          if (g == 0u && l == 0u && m == -1 * static_cast<int>(l))
            debugOutput("    - Adding auxvariables...");

          addAuxVariables(var_name);
        }

        // Add auxkernels.
        if (_current_task == "add_aux_kernel")
        {
          if (g == 0u && l == 0u && m == -1 * static_cast<int>(l))
            debugOutput("    - Adding auxkernels...");

          addAuxKernels(var_name, g, l, m);
        }

        moment_index++;
      }
    }
  }

  if (_current_task == "add_output")
  {
    debugOutput("    - Adding outputs...");

    addOutputs();
  }
}
//------------------------------------------------------------------------------

//------------------------------------------------------------------------------
// Initialize the CGFEM-SAAF scheme.
//------------------------------------------------------------------------------
void
NeutronTransportAction::actSAAFCFEM()
{
  // Loop over all groups.
  for (unsigned int g = 0; g < _num_groups; ++g)
  {
    // Loop over all the flux ordinates.
    for (unsigned int n = 0; n < _num_flux_ordinates; ++n)
    {
      const auto & var_name = _group_angular_fluxes[g][n];

      // Add kernels.
      if (_current_task == "add_kernel")
      {
        if (g == 0u && n == 0u)
          debugOutput("    - Adding kernels...");

        addSAAFKernels(var_name, g, n);

        if (g == _num_groups - 1u && n == _num_flux_ordinates - 1u)
          debugOutput("-----------------------------------------------------",
                      "-----------------------------------------------------");
      }
    }
  }

  // Add Dirac kernels (isotropic point sources). Anisotropic point sources are
  // handled separately.
  if (_current_task == "add_dirac_kernel")
  {
    debugOutput("    - Adding Dirac kernels...");

    addSAAFDiracKernels();
  }
}
//------------------------------------------------------------------------------

//------------------------------------------------------------------------------
// Initialize the DGFEM-upwinding scheme.
//------------------------------------------------------------------------------
void
NeutronTransportAction::actUpwindDFEM()
{
  // Loop over all groups.
  for (unsigned int g = 0; g < _num_groups; ++g)
  {
    // Loop over all the flux ordinates.
    for (unsigned int n = 0; n < _num_flux_ordinates; ++n)
    {
      const auto & var_name = _group_angular_fluxes[g][n];

      // Add kernels.
      if (_current_task == "add_kernel")
      {
        if (g == 0u && n == 0u)
          debugOutput("    - Adding kernels...");

        addDGFEMKernels(var_name, g, n);
      }

      // Add DG kernels.
      if (_current_task == "add_dg_kernel")
      {
        if (g == 0u && n == 0u)
          debugOutput("    - Adding DG kernels...");

        addDGFEMDGKernels(var_name, g, n);

        if (g == _num_groups - 1u && n == _num_flux_ordinates - 1u)
          debugOutput("-----------------------------------------------------",
                      "-----------------------------------------------------");
      }
    }
  }

  // Add Dirac kernels (isotropic point sources). Anisotropic point sources are
  // handled separately.
  if (_current_task == "add_dirac_kernel")
  {
    debugOutput("    - Adding Dirac kernels...");

    addDGFEMDiracKernels();
  }
}
//------------------------------------------------------------------------------

//------------------------------------------------------------------------------
// Functions to add MOOSE objects common to all schemes.
//------------------------------------------------------------------------------
void
NeutronTransportAction::addOutputs()
{
  // Add Exodus
  {
    auto params = _factory.getValidParams("Exodus");

    if (!getParam<bool>("output_angular_fluxes"))
    {
      for (const auto & [g, flux_ordinate] : _group_angular_fluxes)
      {
        for (unsigned int n = 0u; n < flux_ordinate.size(); ++n)
          params.set<std::vector<VariableName>>("hide").emplace_back(flux_ordinate[n]);
      }
    }

    _problem->addOutput("Exodus", "t_out", params);
    debugOutput("      - Adding Output Exodus.");
  } // Exodus
}

void
NeutronTransportAction::addVariable(const std::string & var_name)
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

void
NeutronTransportAction::addBCs(const std::string & var_name, unsigned int g,
                               unsigned int n)
{
  // Add ADSNVacuumBC.
  if (_vacuum_side_sets.size() > 0u)
  {
    auto params = _factory.getValidParams("ADSNVacuumBC");
    params.set<NonlinearVariableName>("variable") = var_name;
    // Ordinate index is required to fetch the neutron direction.
    params.set<unsigned int>("ordinate_index") = n;

    // Apply the parameters for the quadrature rule.
    applyQuadratureParameters(params);

    params.set<std::vector<BoundaryName>>("boundary") = _vacuum_side_sets;

    _problem->addBoundaryCondition("ADSNVacuumBC",
                                   "ADSNVacuumBC_" + var_name,
                                   params);
    debugOutput("      - Adding BC ADSNVacuumBC for the variable "
                + var_name + ".");
  } // ADSNVacuumBC

  // Add ADSNMatSourceBC.
  if (_source_side_sets.size() > 0u)
  {
    auto params = _factory.getValidParams("ADSNMatSourceBC");
    params.set<NonlinearVariableName>("variable") = var_name;
    // Ordinate index is required to fetch the neutron direction.
    params.set<unsigned int>("ordinate_index") = n;

    // Apply the parameters for the quadrature rule.
    applyQuadratureParameters(params);

    params.set<std::vector<BoundaryName>>("boundary") = _source_side_sets;

    _problem->addBoundaryCondition("ADSNMatSourceBC",
                                   "ADSNMatSourceBC" + var_name,
                                   params);
    debugOutput("Adding BC ADSNMatSourceBC for the variable "
                + var_name + ".");
  } // ADSNMatSourceBC

  // Add ADSNReflectiveBC.
  if (_reflective_side_sets.size() > 0u)
  {
    auto params = _factory.getValidParams("ADSNReflectiveBC");
    params.set<NonlinearVariableName>("variable") = var_name;
    // Ordinate index is required to fetch the neutron direction.
    params.set<unsigned int>("ordinate_index") = n;

    // The flux ordinates for this group.
    params.set<std::vector<VariableName>>("psi_ref")
      = _group_angular_fluxes[g];

    // Apply the parameters for the quadrature rule.
    applyQuadratureParameters(params);

    params.set<std::vector<BoundaryName>>("boundary") = _reflective_side_sets;

    _problem->addBoundaryCondition("ADSNReflectiveBC",
                                   "ADSNReflectiveBC" + var_name,
                                   params);
    debugOutput("Adding BC ADSNReflectiveBC for the variable "
                + var_name + ".");
  } // ADSNReflectiveBC
}

// TODO: Initial conditions.
void
NeutronTransportAction::addICs(const std::string & var_name, unsigned int g,
                               unsigned int n)
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
          params.set<std::vector<SubdomainName>>("block")
            = getParam<std::vector<SubdomainName>>("block");
        }

        const auto & const_ic = getParam<std::vector<Real>>("constant_ic");
        if (const_ic.size() == 1)
          params.set<Real>("value") = const_ic[0];
        else if (const_ic.size() == _num_groups)
          params.set<Real>("value") = const_ic[g];
        else
        {
          mooseError("Size of 'constant_ic' does not match the declared number "
                     "of neutron groups.");
        }

        _problem->addInitialCondition("ConstantIC", "ConstantIC_" + var_name, params);
        debugOutput("Adding IC ConstantIC for the variable "
                    + var_name + ".");
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
NeutronTransportAction::addAuxVariables(const std::string & var_name)
{
  auto fe_type = AddVariableAction::feType(_pars);
  auto type = AddVariableAction::variableType(fe_type, false, false);
  auto params = _factory.getValidParams(type);
  params.set<MooseEnum>("order") = fe_type.order.get_order();
  params.set<MooseEnum>("family") = Moose::stringify(fe_type.family);

  if (isParamValid("block"))
  {
    params.set<std::vector<SubdomainName>>("block")
      = getParam<std::vector<SubdomainName>>("block");
  }

  _problem->addAuxVariable(type, var_name, params);
  debugOutput("      - Adding auxvariable " + var_name + ".");
}

void
NeutronTransportAction::addAuxKernels(const std::string & var_name,
                                      unsigned int g,
                                      unsigned int l,
                                      int m)
{
  // Add NeutronFluxMoment.
  {
    InputParameters params = _factory.getValidParams("NeutronFluxMoment");
    params.set<AuxVariableName>("variable") = var_name;
    // Flux moment degree and order.
    params.set<unsigned int>("degree") = l;
    params.set<int>("order") = m;

    // Apply the parameters for the quadrature rule.
    applyQuadratureParameters(params);

    // The flux ordinates for this group.
    params.set<std::vector<VariableName>>("group_flux_ordinates")
      = _group_angular_fluxes[g];

    if (isParamValid("block"))
    {
      params.set<std::vector<SubdomainName>>("block")
        = getParam<std::vector<SubdomainName>>("block");
    }

    _problem->addAuxKernel("NeutronFluxMoment",
                           "NeutronFluxMoment_" + var_name,
                           params);
    debugOutput("      - Adding auxkernel NeutronFluxMoment for the variable "
                + var_name + ".");
  } // NeutronFluxMoment
}
//------------------------------------------------------------------------------

//------------------------------------------------------------------------------
// Functions to add MOOSE objects for the CGFEM-SAAF scheme.
//------------------------------------------------------------------------------
void
NeutronTransportAction::addSAAFKernels(const std::string & var_name,
                                       unsigned int g, unsigned int n)
{
  // Add ADSAAFTimeDerivative.
  if (_exec_type == ExecutionType::Transient)
  {
    auto params = _factory.getValidParams("ADSAAFTimeDerivative");
    params.set<NonlinearVariableName>("variable") = var_name;
    // Group index is required to fetch the group neutron velocity.
    params.set<unsigned int>("group_index") = g;
    // Ordinate index is required to fetch the direction of travel for
    // stabilization.
    params.set<unsigned int>("ordinate_index") = n;

    // Apply the parameters for the quadrature rule.
    applyQuadratureParameters(params);

    if (isParamValid("block"))
    {
      params.set<std::vector<SubdomainName>>("block")
        = getParam<std::vector<SubdomainName>>("block");
    }

    _problem->addKernel("ADSAAFTimeDerivative",
                        "ADSAAFTimeDerivative_" + var_name,
                        params);
    debugOutput("      - Adding kernel ADSAAFTimeDerivative for the variable "
                + var_name + ".");
  } // ADSAAFTimeDerivative

  // Add ADSAAFStreaming.
  {
    auto params = _factory.getValidParams("ADSAAFStreaming");
    params.set<NonlinearVariableName>("variable") = var_name;
    // Group index is required to fetch the group neutron removal cross-section
    // for stabilization.
    params.set<unsigned int>("group_index") = g;
    // Ordinate index is required to fetch the neutron direction.
    params.set<unsigned int>("ordinate_index") = n;

    // Apply the parameters for the quadrature rule.
    applyQuadratureParameters(params);

    if (isParamValid("block"))
    {
      params.set<std::vector<SubdomainName>>("block")
        = getParam<std::vector<SubdomainName>>("block");
    }

    _problem->addKernel("ADSAAFStreaming",
                        "ADSAAFStreaming_" + var_name,
                        params);
    debugOutput("      - Adding kernel ADSAAFStreaming for the variable "
                + var_name + ".");
  } // ADSAAFStreaming

  // Add ADSNRemoval.
  {
    auto params = _factory.getValidParams("ADSNRemoval");
    params.set<NonlinearVariableName>("variable") = var_name;
    // Group index is required to fetch the group neutron removal
    // cross-section.
    params.set<unsigned int>("group_index") = g;

    if (isParamValid("block"))
    {
      params.set<std::vector<SubdomainName>>("block")
        = getParam<std::vector<SubdomainName>>("block");
    }

    _problem->addKernel("ADSNRemoval",
                        "ADSNRemoval_" + var_name,
                        params);
    debugOutput("      - Adding kernel ADSNRemoval for the variable "
                + var_name + ".");
  } // ADSNRemoval

  // Add ADSAAFMaterialSource.
  {
    auto params = _factory.getValidParams("ADSAAFMaterialSource");
    params.set<NonlinearVariableName>("variable") = var_name;
    // Group index and the number of groups are required to fetch the
    // source moments for the proper spectral energy group.
    params.set<unsigned int>("group_index") = g;
    params.set<unsigned int>("num_groups") = _num_groups;
    // Ordinate index is required to fetch the neutron direction.
    params.set<unsigned int>("ordinate_index") = n;

    // Apply the parameters for the quadrature rule.
    applyQuadratureParameters(params);

    if (isParamValid("block"))
    {
      params.set<std::vector<SubdomainName>>("block")
        = getParam<std::vector<SubdomainName>>("block");
    }

    _problem->addKernel("ADSAAFMaterialSource",
                        "ADSAAFMaterialSource_" + var_name,
                        params);
    debugOutput("      - Adding kernel ADSAAFMaterialSource for the variable "
                + var_name + ".");
  } // ADSAAFMaterialSource

  // Only add scattering kernels if debug doesn't disable them.
  if (!getParam<bool>("debug_disable_scattering"))
  {
    // Debug option to disable the source iteration solver.
    if (!getParam<bool>("debug_disable_source_iteration"))
    {
      // Compute the scattering evaluation with source iteration.
      mooseError("Scattering iteration is not supported and is a work in "
                 "progress.");

      // Add ADSAAFExternalScattering. This kernel is only enabled if there's
      // more than one group.
      if (_num_groups > 1u)
      {
        auto params = _factory.getValidParams("ADSAAFExternalScattering");
        params.set<NonlinearVariableName>("variable") = var_name;
        // Group index and the number of groups are required to fetch the
        // source moments for the proper spectral energy group.
        params.set<unsigned int>("group_index") = g;
        params.set<unsigned int>("num_groups") = _num_groups;
        // Ordinate index is required to fetch the neutron direction.
        params.set<unsigned int>("ordinate_index") = n;

        // Apply the parameters for the quadrature rule.
        applyQuadratureParameters(params);

        // Copy all of the group flux moment names into the moment variable
        // parameter.
        auto & moment_names = params.set<std::vector<VariableName>>("group_flux_moments");
        for (unsigned int g_prime = 0; g_prime < _num_groups; ++g_prime)
        {
          std::copy(_group_flux_moments[g_prime].begin(),
                    _group_flux_moments[g_prime].end(),
                    std::back_inserter(moment_names));
        }

        if (isParamValid("block"))
        {
          params.set<std::vector<SubdomainName>>("block")
            = getParam<std::vector<SubdomainName>>("block");
        }

        _problem->addKernel("ADSAAFExternalScattering",
                            "ADSAAFExternalScattering_" + var_name,
                            params);
        debugOutput("      - Adding kernel ADSAAFExternalScattering for the variable "
                    + var_name + ".");
      } // ADSAAFExternalScattering

      // Add ADSAAFInternalScattering.
      {
        auto params = _factory.getValidParams("ADSAAFInternalScattering");
        params.set<NonlinearVariableName>("variable") = var_name;
        // Group index and the number of groups are required to fetch the
        // source moments for the proper spectral energy group.
        params.set<unsigned int>("group_index") = g;
        params.set<unsigned int>("num_groups") = _num_groups;
        // Ordinate index is required to fetch the neutron direction.
        params.set<unsigned int>("ordinate_index") = n;

        // Apply the parameters for the quadrature rule.
        applyQuadratureParameters(params);

        // Copy all of the in-group flux moment names into the moment
        // variable parameter.
        params.set<std::vector<VariableName>>("within_group_flux_moments")
          = _group_flux_moments[g];

        if (isParamValid("block"))
        {
          params.set<std::vector<SubdomainName>>("block")
            = getParam<std::vector<SubdomainName>>("block");
        }

        _problem->addKernel("ADSAAFInternalScattering",
                            "ADSAAFInternalScattering_" + var_name,
                            params);
        debugOutput("      - Adding kernel ADSAAFInternalScattering for the variable "
                    + var_name + ".");
      } // ADSAAFInternalScattering
    }
    else
    {
      // Computes the scattering evaluation without source iteration.
      // Add ADSAAFScattering.
      {
        auto params = _factory.getValidParams("ADSAAFScattering");
        params.set<NonlinearVariableName>("variable") = var_name;
        // Group index and the number of groups are required to fetch the
        // source moments for the proper spectral energy group.
        params.set<unsigned int>("group_index") = g;
        params.set<unsigned int>("num_groups") = _num_groups;
        // Ordinate index is required to fetch the neutron direction.
        params.set<unsigned int>("ordinate_index") = n;
        // Maximum scattering anisotropy.
        params.set<unsigned int>("max_anisotropy") = _max_eval_anisotropy;

        // Apply the parameters for the quadrature rule.
        applyQuadratureParameters(params);

        // Copy all of the group flux ordinate names into the variable
        // parameter.
        auto & ordinate_names = params.set<std::vector<VariableName>>("group_flux_ordinates");
        for (unsigned int g_prime = 0; g_prime < _num_groups; ++g_prime)
        {
          std::copy(_group_angular_fluxes[g_prime].begin(),
                    _group_angular_fluxes[g_prime].end(),
                    std::back_inserter(ordinate_names));
        }

        if (isParamValid("block"))
        {
          params.set<std::vector<SubdomainName>>("block")
            = getParam<std::vector<SubdomainName>>("block");
        }

        _problem->addKernel("ADSAAFScattering",
                            "ADSAAFScattering_" + var_name,
                            params);
        debugOutput("      - Adding kernel ADSAAFScattering for the variable "
                    + var_name + ".");
      } // ADSAAFScattering
    }
  }
}

void
NeutronTransportAction::addSAAFDiracKernels()
{
  // Add DFEMIsoPointSource.
  {
    const auto & source_locations
      = getParam<std::vector<Point>>("point_source_locations");
    const auto & source_intensities
      = getParam<std::vector<Real>>("point_source_intensities");
    const auto & source_groups
      = getParam<std::vector<unsigned int>>("point_source_groups");

    if (source_locations.size() != source_intensities.size() ||
        source_intensities.size() != source_groups.size())
    {
      mooseWarning("The number of provided parameters for the isotropic point "
                   "sources do not match. Some sources will be ignored and "
                   "others may be mismatched between their locations, "
                   "intensities, and emission groups.");
    }

    // Loop over the vector of group indices to bin the source indices by the
    // group they emit into.
    std::unordered_map<unsigned int, std::vector<unsigned int>> group_map;
    const unsigned int num_sources = std::min(std::min(source_locations.size(),
                                                       source_intensities.size()),
                                                       source_groups.size());
    for (unsigned int i = 0; i < num_sources; ++i)
    {
      if (source_groups[i] > _num_groups)
      {
        mooseWarning("Group " + Moose::stringify(source_groups[i]) + " exceeds "
                     "the number of groups requested. Ignoring this source.");
        continue;
      }

      if (group_map.count(source_groups[i]) > 0u)
        group_map.emplace(source_groups[i] - 1u, std::vector<unsigned int>());
      else
        group_map[source_groups[i] - 1u].emplace_back(i);
    }

    // Loop over the sorted groups and assign dirac kernels to each flux group.
    for (const auto & [g, mapping] : group_map)
    {
      for (unsigned int n = 0; n < _num_flux_ordinates; ++n)
      {
        const auto & var_name = _group_angular_fluxes[g][n];

        auto params = _factory.getValidParams("SAAFIsoPointSource");
        params.set<NonlinearVariableName>("variable") = var_name;
        params.set<MooseEnum>("point_not_found_behavior")
          = MooseEnum("ERROR WARNING IGNORE", "WARNING");
        // Group index is required to fetch the group neutron removal cross-section
        // for stabilization.
        params.set<unsigned int>("group_index") = g;
        // Ordinate index is required to fetch the neutron direction.
        params.set<unsigned int>("ordinate_index") = n;

        // Apply the parameters for the quadrature rule.
        applyQuadratureParameters(params);

        // Set the intensities and points for the group g.
        for (const auto & index : mapping)
        {
          params.set<std::vector<Point>>("points").emplace_back(source_locations[index]);
          params.set<std::vector<Real>>("intensities").emplace_back(source_intensities[index]);
        }

        if (isParamValid("block"))
        {
          params.set<std::vector<SubdomainName>>("block")
            = getParam<std::vector<SubdomainName>>("block");
        }

        _problem->addDiracKernel("SAAFIsoPointSource",
                                 "SAAFIsoPointSource_" + var_name,
                                 params);
        debugOutput("      - Adding Dirac kernel SAAFIsoPointSource for the "
                    "variable " + var_name + ".");
      }
    }
  } // SAAFIsoPointSource
}
//------------------------------------------------------------------------------

//------------------------------------------------------------------------------
// Functions to add MOOSE objects for the DGFEM-upwinding scheme.
//------------------------------------------------------------------------------
void
NeutronTransportAction::addDGFEMKernels(const std::string & var_name, unsigned int g,
                                        unsigned int n)
{
  // Add ADDFEMTimeDerivative.
  if (_exec_type == ExecutionType::Transient)
  {
    auto params = _factory.getValidParams("ADDFEMTimeDerivative");
    params.set<NonlinearVariableName>("variable") = var_name;
    // Group index is required to fetch the group neutron velocity.
    params.set<unsigned int>("group_index") = g;

    if (isParamValid("block"))
    {
      params.set<std::vector<SubdomainName>>("block")
        = getParam<std::vector<SubdomainName>>("block");
    }

    _problem->addKernel("ADDFEMTimeDerivative",
                        "ADDFEMTimeDerivative_" + var_name,
                        params);
    debugOutput("      - Adding kernel ADDFEMTimeDerivative for the variable "
                + var_name + ".");
  } // ADDFEMTimeDerivative

  // Add ADDFEMStreaming.
  {
    auto params = _factory.getValidParams("ADDFEMStreaming");
    params.set<NonlinearVariableName>("variable") = var_name;
    // Ordinate index is required to fetch the neutron direction.
    params.set<unsigned int>("ordinate_index") = n;

    // Apply the parameters for the quadrature rule.
    applyQuadratureParameters(params);

    if (isParamValid("block"))
    {
      params.set<std::vector<SubdomainName>>("block")
        = getParam<std::vector<SubdomainName>>("block");
    }

    _problem->addKernel("ADDFEMStreaming",
                        "ADDFEMStreaming_" + var_name,
                        params);
    debugOutput("      - Adding kernel ADDFEMStreaming for the variable "
                + var_name + ".");
  } // ADDFEMStreaming

  // Add ADSNRemoval.
  {
    auto params = _factory.getValidParams("ADSNRemoval");
    params.set<NonlinearVariableName>("variable") = var_name;
    // Group index is required to fetch the group neutron removal
    // cross-section.
    params.set<unsigned int>("group_index") = g;

    if (isParamValid("block"))
    {
      params.set<std::vector<SubdomainName>>("block")
        = getParam<std::vector<SubdomainName>>("block");
    }

    _problem->addKernel("ADSNRemoval",
                        "ADSNRemoval_" + var_name,
                        params);
    debugOutput("      - Adding kernel ADSNRemoval for the variable "
                + var_name + ".");
  } // ADSNRemoval

  // Add ADDFEMMaterialSource.
  {
    auto params = _factory.getValidParams("ADDFEMMaterialSource");
    params.set<NonlinearVariableName>("variable") = var_name;
    // Group index and the number of groups are required to fetch the
    // source moments for the proper spectral energy group.
    params.set<unsigned int>("group_index") = g;
    params.set<unsigned int>("num_groups") = _num_groups;
    // Ordinate index is required to fetch the neutron direction.
    params.set<unsigned int>("ordinate_index") = n;

    // Apply the parameters for the quadrature rule.
    applyQuadratureParameters(params);

    if (isParamValid("block"))
    {
      params.set<std::vector<SubdomainName>>("block")
        = getParam<std::vector<SubdomainName>>("block");
    }

    _problem->addKernel("ADDFEMMaterialSource",
                        "ADDFEMMaterialSource_" + var_name,
                        params);
    debugOutput("      - Adding kernel ADDFEMMaterialSource for the variable "
                + var_name + ".");
  } // ADDFEMMaterialSource

  // Only add scattering kernels if debug doesn't disable them.
  if (!getParam<bool>("debug_disable_scattering"))
  {
    // Debug option to disable the source iteration solver.
    if (!getParam<bool>("debug_disable_source_iteration"))
    {
      // Compute the scattering evaluation with source iteration.
      mooseError("Scattering iteration is not supported and is a work in "
                 "progress.");

      // Add ADDFEMExternalScattering. This kernel is only enabled if there's
      // more than one group.
      if (_num_groups > 1u)
      {
        auto params = _factory.getValidParams("ADDFEMExternalScattering");
        params.set<NonlinearVariableName>("variable") = var_name;
        // Group index and the number of groups are required to fetch the
        // source moments for the proper spectral energy group.
        params.set<unsigned int>("group_index") = g;
        params.set<unsigned int>("num_groups") = _num_groups;
        // Ordinate index is required to fetch the neutron direction.
        params.set<unsigned int>("ordinate_index") = n;

        // Apply the parameters for the quadrature rule.
        applyQuadratureParameters(params);

        // Copy all of the group flux moment names into the moment variable
        // parameter.
        auto & moment_names = params.set<std::vector<VariableName>>("group_flux_moments");
        for (unsigned int g_prime = 0; g_prime < _num_groups; ++g_prime)
        {
          std::copy(_group_flux_moments[g_prime].begin(),
                    _group_flux_moments[g_prime].end(),
                    std::back_inserter(moment_names));
        }

        if (isParamValid("block"))
        {
          params.set<std::vector<SubdomainName>>("block")
            = getParam<std::vector<SubdomainName>>("block");
        }

        _problem->addKernel("ADDFEMExternalScattering",
                            "ADDFEMExternalScattering_" + var_name,
                            params);
        debugOutput("      - Adding kernel ADDFEMExternalScattering for the variable "
                    + var_name + ".");
      } // ADDFEMExternalScattering

      // Add ADDFEMInternalScattering.
      {
        auto params = _factory.getValidParams("ADDFEMInternalScattering");
        params.set<NonlinearVariableName>("variable") = var_name;
        // Group index and the number of groups are required to fetch the
        // source moments for the proper spectral energy group.
        params.set<unsigned int>("group_index") = g;
        params.set<unsigned int>("num_groups") = _num_groups;
        // Ordinate index is required to fetch the neutron direction.
        params.set<unsigned int>("ordinate_index") = n;

        // Apply the parameters for the quadrature rule.
        applyQuadratureParameters(params);

        // Copy all of the in-group flux moment names into the moment
        // variable parameter.
        params.set<std::vector<VariableName>>("within_group_flux_moments")
          = _group_flux_moments[g];

        if (isParamValid("block"))
        {
          params.set<std::vector<SubdomainName>>("block")
            = getParam<std::vector<SubdomainName>>("block");
        }

        _problem->addKernel("ADDFEMInternalScattering",
                            "ADDFEMInternalScattering_" + var_name,
                            params);
        debugOutput("      - Adding kernel ADDFEMInternalScattering for the variable "
                    + var_name + ".");
      } // ADDFEMInternalScattering
    }
    else
    {
      // Computes the scattering evaluation without source iteration.
      // Add ADDFEMScattering.
      {
        auto params = _factory.getValidParams("ADDFEMScattering");
        params.set<NonlinearVariableName>("variable") = var_name;
        // Group index and the number of groups are required to fetch the
        // source moments for the proper spectral energy group.
        params.set<unsigned int>("group_index") = g;
        params.set<unsigned int>("num_groups") = _num_groups;
        // Ordinate index is required to fetch the neutron direction.
        params.set<unsigned int>("ordinate_index") = n;
        // Maximum scattering anisotropy.
        params.set<unsigned int>("max_anisotropy") = _max_eval_anisotropy;

        // Apply the parameters for the quadrature rule.
        applyQuadratureParameters(params);

        // Copy all of the group flux ordinate names into the variable
        // parameter.
        auto & ordinate_names = params.set<std::vector<VariableName>>("group_flux_ordinates");
        for (unsigned int g_prime = 0; g_prime < _num_groups; ++g_prime)
        {
          std::copy(_group_angular_fluxes[g_prime].begin(),
                    _group_angular_fluxes[g_prime].end(),
                    std::back_inserter(ordinate_names));
        }

        if (isParamValid("block"))
        {
          params.set<std::vector<SubdomainName>>("block")
            = getParam<std::vector<SubdomainName>>("block");
        }

        _problem->addKernel("ADDFEMScattering",
                            "ADDFEMScattering_" + var_name,
                            params);
        debugOutput("      - Adding kernel ADDFEMScattering for the variable "
                    + var_name + ".");
      } // ADDFEMScattering
    }
  }
}

void
NeutronTransportAction::addDGFEMDGKernels(const std::string & var_name,
                                          unsigned int g, unsigned int n)
{
  // Add ADDFEMUpwinding.
  {
    auto params = _factory.getValidParams("ADDFEMUpwinding");
    params.set<NonlinearVariableName>("variable") = var_name;
    // Ordinate index is required to fetch the neutron direction.
    params.set<unsigned int>("ordinate_index") = n;

    // Apply the parameters for the quadrature rule.
    applyQuadratureParameters(params);

    if (isParamValid("block"))
    {
      params.set<std::vector<SubdomainName>>("block")
        = getParam<std::vector<SubdomainName>>("block");
    }

    _problem->addDGKernel("ADDFEMUpwinding",
                          "ADDFEMUpwinding_" + var_name,
                          params);
    debugOutput("      - Adding DG kernel ADDFEMUpwinding for the variable "
                + var_name + ".");
  } // ADDFEMUpwinding
}

void
NeutronTransportAction::addDGFEMDiracKernels()
{
  // Add DFEMIsoPointSource.
  {
    const auto & source_locations
      = getParam<std::vector<Point>>("point_source_locations");
    const auto & source_intensities
      = getParam<std::vector<Real>>("point_source_intensities");
    const auto & source_groups
      = getParam<std::vector<unsigned int>>("point_source_groups");

    if (source_locations.size() != source_intensities.size() ||
        source_intensities.size() != source_groups.size())
    {
      mooseWarning("The number of provided parameters for the isotropic point "
                   "sources do not match. Some sources will be ignored and "
                   "others may be mismatched between their locations, "
                   "intensities, and emission groups.");
    }

    // Loop over the vector of group indices to bin the source indices by the
    // group they emit into.
    std::unordered_map<unsigned int, std::vector<unsigned int>> group_map;
    const unsigned int num_sources = std::min(std::min(source_locations.size(),
                                                       source_intensities.size()),
                                                       source_groups.size());
    for (unsigned int i = 0; i < num_sources; ++i)
    {
      if (source_groups[i] > _num_groups)
      {
        mooseWarning("Group " + Moose::stringify(source_groups[i]) + " exceeds "
                     "the number of groups requested. Ignoring this source.");
        continue;
      }

      if (group_map.count(source_groups[i]) > 0u)
        group_map.emplace(source_groups[i] - 1u, std::vector<unsigned int>());
      else
        group_map[source_groups[i] - 1u].emplace_back(i);
    }

    // Loop over the sorted groups and assign dirac kernels to each flux group.
    for (const auto & [g, mapping] : group_map)
    {
      for (unsigned int n = 0; n < _num_flux_ordinates; ++n)
      {
        const auto & var_name = _group_angular_fluxes[g][n];

        auto params = _factory.getValidParams("DFEMIsoPointSource");
        params.set<NonlinearVariableName>("variable") = var_name;
        params.set<MooseEnum>("point_not_found_behavior")
          = MooseEnum("ERROR WARNING IGNORE", "WARNING");

        // Apply the parameters for the quadrature rule.
        applyQuadratureParameters(params);

        // Set the intensities and points for the group g.
        for (const auto & index : mapping)
        {
          params.set<std::vector<Point>>("points").emplace_back(source_locations[index]);
          params.set<std::vector<Real>>("intensities").emplace_back(source_intensities[index]);
        }

        if (isParamValid("block"))
        {
          params.set<std::vector<SubdomainName>>("block")
            = getParam<std::vector<SubdomainName>>("block");
        }

        _problem->addDiracKernel("DFEMIsoPointSource",
                                 "DFEMIsoPointSource_" + var_name,
                                 params);
        debugOutput("      - Adding Dirac kernel DFEMIsoPointSource for the "
                    "variable " + var_name + ".");
      }
    }
  } // DFEMIsoPointSource
}
//------------------------------------------------------------------------------
