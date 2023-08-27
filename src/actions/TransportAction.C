#include "TransportAction.h"

#include "Factory.h"
#include "Parser.h"
#include "NonlinearSystemBase.h"
#include "AuxiliarySystem.h"
#include "FEProblemBase.h"
#include "ConsoleStream.h"

#include "Conversion.h"
#include "MooseTypes.h"
#include "FEProblem.h"

#include "AddVariableAction.h"
#include "AddOutputAction.h"
#include "ActionWarehouse.h"

#include "libmesh/string_to_enum.h"
#include "libmesh/fe_type.h"

#include "RayKernelBase.h"

// All schemes.
registerMooseAction("GnatApp", TransportAction, "add_variable");
registerMooseAction("GnatApp", TransportAction, "add_kernel");
registerMooseAction("GnatApp", TransportAction, "add_dirac_kernel");
registerMooseAction("GnatApp", TransportAction, "add_bc");
registerMooseAction("GnatApp", TransportAction, "add_ic");
registerMooseAction("GnatApp", TransportAction, "add_aux_variable");
registerMooseAction("GnatApp", TransportAction, "add_aux_kernel");
registerMooseAction("GnatApp", TransportAction, "add_user_object");

// Picard iteration.
registerMooseAction("GnatApp", TransportAction, "add_transfer");
// Required for conservative transfers between MOOSE applications.
registerMooseAction("GnatApp", TransportAction, "add_postprocessor");

// Restart.
registerMooseAction("GnatApp", TransportAction, "check_copy_nodal_vars");
registerMooseAction("GnatApp", TransportAction, "copy_nodal_vars");

InputParameters
TransportAction::validParams()
{
  auto params = GnatBaseAction::validParams();
  params.addClassDescription("This action adds all of the required variables, "
                             "kernels, boundary conditions, initial conditions "
                             "and auxiliary systems required to solve "
                             "source-driven multi-group neutral particle transport "
                             "problems with Gnat.");

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
  params.addParam<bool>("output_angular_fluxes",
                        false,
                        "Whether the angular flux ordinates should be written "
                        "to the exodus file or not.");
  params.addParamNamesToGroup("scaling angular_flux_names flux_moment_names output_angular_fluxes ",
                              "Variable");

  //----------------------------------------------------------------------------
  // Basic parameters for the transport simulation.
  params.addRequiredParam<MooseEnum>("scheme",
                                     MooseEnum("saaf_cfem diffusion_cfem flux_moment_transfer"),
                                     "The discretization and stabilization "
                                     "scheme for the transport equation.");
  params.addRequiredParam<MooseEnum>(
      "particle_type",
      MooseEnum("neutron photon"),
      "The type of particle to be consuming material property data.");
  params.addRequiredParam<unsigned int>("num_groups",
                                        "The number of spectral energy groups in the "
                                        "problem.");
  params.addParam<unsigned int>("max_anisotropy",
                                0,
                                "The maximum degree of anisotropy to evaluate. "
                                "Defaults to 0 for isotropic scattering.");
  params.addParam<bool>(
      "eigen",
      false,
      "Whether or not the transport simulation should prepare for an eigenvalue calculation. This "
      "parameters is only valid for steady-state simulations with the particle type set to "
      "'neutron'. This parameter is ignored otherwise.");
  params.addParam<bool>(
      "init_from_file", false, "If the simulation should be initialized from a file or not.");
  params.addParam<bool>(
      "use_scattering_jacobians", false, "Whether or not to use hand-coded scattering Jacobians.");

  params.addParamNamesToGroup("eigen max_anisotropy block use_scattering_jacobians init_from_file",
                              "Simulation");

  //----------------------------------------------------------------------------
  // Quadrature parameters.
  params.addRangeCheckedParam<unsigned int>("n_polar",
                                            3,
                                            "n_polar > 0",
                                            "Number of Legendre polar "
                                            "quadrature points in a single "
                                            "octant of the unit sphere. "
                                            "Defaults to 3.");
  params.addRangeCheckedParam<unsigned int>("n_azimuthal",
                                            3,
                                            "n_azimuthal > 0",
                                            "Number of Chebyshev azimuthal "
                                            "quadrature points in a single "
                                            "octant of the unit sphere. "
                                            "Defaults to 3.");
  params.addParam<MooseEnum>("major_axis",
                             MooseEnum("x y z", "x"),
                             "Major axis of the angular quadrature. Allows the "
                             "polar angular quadrature to align with a cartesian "
                             "axis with minimal heterogeneity. Default is the "
                             "x-axis. This parameter is only applied in 3D "
                             "cartesian problems.");
  params.addParamNamesToGroup("n_polar n_azimuthal major_axis", "Quadrature");

  //----------------------------------------------------------------------------
  // Source-driven problem parameters.
  params.addParam<bool>("scale_sources",
                        false,
                        "If the external fixed sources should be scaled such that each source "
                        "moment is divided by the maximum source moment.");
  params.addParamNamesToGroup("scale_sources", "Fixed Source");

  //----------------------------------------------------------------------------
  // Boundary conditions.
  params.addParam<std::vector<BoundaryName>>("vacuum_boundaries",
                                             "The boundaries to apply vacuum "
                                             "boundary conditions.");
  params.addParam<std::vector<BoundaryName>>("source_boundaries",
                                             "The boundaries to apply incoming "
                                             "flux boundary conditions.");
  params.addParam<std::vector<BoundaryName>>(
      "current_boundaries",
      "The boundaries to apply the current boundary conditions to. This is a specialization of the "
      "source boundary condition for a particle direction equal to the surface normal.");
  params.addParam<std::vector<BoundaryName>>("reflective_boundaries",
                                             "The boundaries to apply reflective "
                                             "boundary conditions.");

  params.addParam<std::vector<std::vector<Real>>>(
      "boundary_source_moments",
      "A double vector containing the external source moments for "
      "all boundaries. The exterior vector must correspond with the surface source boundary "
      "conditions provided in 'source_boundaries'.");
  params.addParam<std::vector<unsigned int>>(
      "boundary_source_anisotropy",
      "The degree of anisotropy of the boundary source moments. The exterior vector must "
      "correspond with the surface source boundary "
      "conditions provided in 'source_boundaries'.");

  params.addParam<std::vector<std::vector<Real>>>(
      "boundary_currents",
      "A double vector containing the external currents for all boundaries. The exterior vector "
      "must correspond with the surface current boundary conditions provided in "
      "'current_boundaries'.");
  params.addParam<std::vector<unsigned int>>(
      "boundary_current_anisotropy",
      "The degree of anisotropy to be applied to the boundary currents. The vector must correspond "
      "with the surface current boundary conditions provided in 'current_boundaries'.");

  params.addParamNamesToGroup(
      "vacuum_boundaries source_boundaries current_boundaries reflective_boundaries "
      "boundary_source_moments boundary_source_anisotropy boundary_currents "
      "boundary_current_anisotropy",
      "Boundary Condition");

  //----------------------------------------------------------------------------
  // Volumetric sources.
  params.addParam<std::vector<SubdomainName>>("volumetric_source_blocks",
                                              "The list of blocks (ids or "
                                              "names) that host a volumetric source.");
  params.addParam<std::vector<std::vector<Real>>>(
      "volumetric_source_moments",
      "A double vector containing a list of external source moments for all volumetric particle "
      "sources. The external vector should correspond with the order of "
      "'volumetric_source_blocks'.");
  params.addParam<std::vector<unsigned int>>(
      "volumetric_source_anisotropies",
      "The anisotropies of the volumetric sources. The vector should correspond with the order of "
      "'volumetric_source_blocks'");
  params.addParamNamesToGroup(
      "volumetric_source_blocks volumetric_source_moments volumetric_source_anisotropies",
      "Volumetric Source");

  //----------------------------------------------------------------------------
  // Point sources.
  params.addParam<std::vector<Point>>("point_source_locations",
                                      "The locations of all isotropic "
                                      "point sources in the problem "
                                      "space.");
  params.addParam<std::vector<std::vector<Real>>>(
      "point_source_moments",
      "A double vector containing a list of external source moments for all point particle "
      "sources. The external vector should correspond with the order of "
      "'point_source_locations'.");
  params.addParam<std::vector<unsigned int>>(
      "point_source_anisotropies",
      "The anisotropies of the point sources. The vector should correspond with the order of "
      "'point_source_locations'");
  params.addParamNamesToGroup("point_source_locations point_source_moments "
                              "point_source_anisotropies",
                              "Point Source");

  //----------------------------------------------------------------------------
  // Initial conditions.
  params.addParam<MooseEnum>("ic_type",
                             MooseEnum("constant function file", "constant"),
                             "The type of initial condition to use (if "
                             "multiple are provided). Defaults to constant "
                             "initial conditions.");
  params.addParam<std::vector<Real>>("constant_ic",
                                     std::vector<Real>(0.0),
                                     "A constant initial condition for the "
                                     "angular fluxes.");
  params.addParamNamesToGroup("ic_type constant_ic", "Initial Condition");
  // TODO: Init from function.

  //----------------------------------------------------------------------------
  // Parameters for multi-app (Picard) coupling.
  params.addParam<MultiAppName>(
      "from_multi_app", "", "The name of the multi-app to pull the source flux moments.");
  params.addParam<std::string>(
      "from_flux_moment_names", "flux_moment", "The names of the source flux moments.");
  params.addParam<std::vector<SubdomainName>>("from_blocks",
                                              "The list of blocks (ids or "
                                              "names) that we are pulling from.");
  params.addParam<bool>(
      "use_copy",
      false,
      "Whether a MultiAppCopyTransfer should be used for the transfer scheme or not.");
  params.addParam<bool>("transfer_to_fv",
                        false,
                        "Whether the transfer variable should be a finite volume variable or not. "
                        "Useful for coupling to finite volume fields.");
  params.addParam<bool>(
      "use_conservative_transfers",
      false,
      "Whether this transport action should pull flux moments using conservative transfers.");
  params.addParam<bool>("is_conservative_transfer_src",
                        false,
                        "Whether this transport action is providing flux moments to a sub/parent "
                        "app using conservative transfers. Setting this option to 'true' adds "
                        "post-processors to ensure flux moments are conservative.");
  params.addParamNamesToGroup(
      "from_multi_app from_flux_moment_names from_blocks use_copy "
      "transfer_to_fv use_conservative_transfers is_conservative_transfer_src",
      "MultiApp");

  //----------------------------------------------------------------------------
  // Parameters for a multi-app provided uncollided flux treatment.
  params.addParam<MultiAppName>("uncollided_from_multi_app",
                                "",
                                "The name of the multi-app to pull the uncollided flux moments.");
  params.addParam<std::string>("from_uncollided_flux_moment_names",
                               "uncollided_flux_moment",
                               "The names of the source flux moments.");
  params.addParam<std::vector<SubdomainName>>("uncollided_from_blocks",
                                              "The list of blocks (ids or "
                                              "names) that we are pulling from.");

  //----------------------------------------------------------------------------
  // Parameters for debugging.
  params.addParam<MooseEnum>("debug_verbosity",
                             MooseEnum("level0 level1", "level1"),
                             "How verbose the debug output of the transport "
                             "system should be. level0 is fully verbose. "
                             "level1 outputs less debugging information.");
  params.addParam<bool>(
      "debug_disable_scattering", false, "Debug option to disable scattering evaluation.");
  params.addParam<bool>(
      "debug_disable_fission", true, "Debug option to disable fission evaluation.");
  params.addParam<bool>(
      "debug_disable_source_iteration", true, "Debug option to disable source iteration.");
  params.addParamNamesToGroup("debug_verbosity debug_disable_scattering debug_disable_fission "
                              "debug_disable_source_iteration",
                              "Debugging");

  params.addParamNamesToGroup("family order num_groups scheme particle_type", "Required");

  return params;
}

TransportAction::TransportAction(const InputParameters & params)
  : GnatBaseAction(params),
    _transport_scheme(getParam<MooseEnum>("scheme").getEnum<TransportScheme>()),
    _particle(getParam<MooseEnum>("particle_type").getEnum<Particletype>()),
    _is_eigen(getParam<bool>("eigen")),
    _n_l(getParam<unsigned int>("n_polar")),
    _n_c(getParam<unsigned int>("n_azimuthal")),
    _num_flux_ordinates(0u),
    _num_group_moments(0u),
    _num_groups(getParam<unsigned int>("num_groups")),
    _max_eval_anisotropy(getParam<unsigned int>("max_anisotropy")),
    _flux_moment_name(getParam<std::string>("flux_moment_names")),
    _angular_flux_name(getParam<std::string>("angular_flux_names")),
    _vacuum_side_sets(getParam<std::vector<BoundaryName>>("vacuum_boundaries")),
    _source_side_sets(getParam<std::vector<BoundaryName>>("source_boundaries")),
    _current_side_sets(getParam<std::vector<BoundaryName>>("current_boundaries")),
    _reflective_side_sets(getParam<std::vector<BoundaryName>>("reflective_boundaries")),
    _point_source_locations(getParam<std::vector<Point>>("point_source_locations")),
    _point_source_moments(getParam<std::vector<std::vector<Real>>>("point_source_moments")),
    _point_source_anisotropy(getParam<std::vector<unsigned int>>("point_source_anisotropies")),
    _boundary_source_moments(getParam<std::vector<std::vector<Real>>>("boundary_source_moments")),
    _boundary_source_anisotropy(getParam<std::vector<unsigned int>>("boundary_source_anisotropy")),
    _boundary_currents(getParam<std::vector<std::vector<Real>>>("boundary_currents")),
    _boundary_current_anisotropy(
        getParam<std::vector<unsigned int>>("boundary_current_anisotropy")),
    _volumetric_source_blocks(getParam<std::vector<SubdomainName>>("volumetric_source_blocks")),
    _volumetric_source_moments(
        getParam<std::vector<std::vector<Real>>>("volumetric_source_moments")),
    _volumetric_source_anisotropy(
        getParam<std::vector<unsigned int>>("volumetric_source_anisotropies")),
    _from_multi_app_name(getParam<MultiAppName>("from_multi_app")),
    _from_subdomain_ids(getParam<std::vector<SubdomainName>>("from_blocks")),
    _source_flux_moment_names(getParam<std::string>("from_flux_moment_names")),
    _uncollided_from_multi_app_name(getParam<MultiAppName>("uncollided_from_multi_app")),
    _uncollided_from_subdomain_ids(getParam<std::vector<SubdomainName>>("uncollided_from_blocks")),
    _uncollided_source_flux_moment_names(
        getParam<std::string>("from_uncollided_flux_moment_names")),
    _using_uncollided(_uncollided_from_multi_app_name != ""),
    _source_scale_factor(0.0),
    _var_init(false)
{
  if (getParam<bool>("scale_sources"))
  {
    // Find the maximum source moment.
    for (const auto & point_moments : _point_source_moments)
      for (auto & point_moment : point_moments)
        _source_scale_factor = std::max(_source_scale_factor, point_moment);

    for (const auto & surf_moments : _boundary_source_moments)
      for (const auto & surf_moment : surf_moments)
        _source_scale_factor = std::max(_source_scale_factor, surf_moment);

    for (const auto & surf_currents : _boundary_currents)
      for (const auto & surf_current : surf_currents)
        _source_scale_factor = std::max(_source_scale_factor, surf_current);

    for (const auto & volume_moments : _volumetric_source_moments)
      for (const auto & volume_moment : volume_moments)
        _source_scale_factor = std::max(_source_scale_factor, volume_moment);

    // Divide all source moments by the maximum. This will be reverted when the scaled flux moments
    // are calculated.
    for (auto & point_moments : _point_source_moments)
      for (auto & point_moment : point_moments)
        point_moment /= _source_scale_factor;

    for (auto & surf_moments : _boundary_source_moments)
      for (auto & surf_moment : surf_moments)
        surf_moment /= _source_scale_factor;

    for (auto & surf_currents : _boundary_currents)
      for (auto & surf_current : surf_currents)
        surf_current /= _source_scale_factor;

    for (auto & volume_moments : _volumetric_source_moments)
      for (auto & volume_moment : volume_moments)
        volume_moment /= _source_scale_factor;
  }

  if (_using_uncollided && _transport_scheme != TransportScheme::SAAFCFEM)
    mooseWarning("Uncollided flux corrections only work for discrete ordinates transport schemes. "
                 "The uncollided flux moments will not be used.");

  if (_using_uncollided && _max_eval_anisotropy > 0u)
    mooseWarning("Currently uncollided flux treatments only support scalar flux moments. Higher "
                 "order / degree momement support is planned in the future.");
}

void
TransportAction::act()
{
  // Initialize common members.
  actCommon();

  // Act functions for different schemes.
  switch (_transport_scheme)
  {
    case TransportScheme::SAAFCFEM:
      actSAAFCFEM();

      if (_using_uncollided && _var_init)
        actUncollided();
      break;
    case TransportScheme::DiffusionApprox:
      actDiffusion();
      break;
    case TransportScheme::FluxMomentTransfer:
      actTransfer();
      break;
    default:
      break;
  }
}

//------------------------------------------------------------------------------
// Function to initialize SN quadrature parameters.
//------------------------------------------------------------------------------
void
TransportAction::applyQuadratureParameters(InputParameters & params)
{
  params.set<UserObjectName>("aq") = "AQProvider_" + name();
}
//------------------------------------------------------------------------------

//------------------------------------------------------------------------------
// Function to initialize common parameters for all schemes.
//------------------------------------------------------------------------------
void
TransportAction::actCommon()
{
  // Setup for all actions the TransportAction action performs.
  // Initialize base action parameters.
  if (_problem && !_var_init)
    initializeBase();

  // Fission warnings.
  if (_problem && !_var_init)
  {
    if (!getParam<bool>("debug_disable_fission") && _particle == Particletype::Neutron)
    {
      std::string message(
          "Gnat is not designed for the analysis of fissile systems. Calculations including "
          "neutron multiplication have been provided for the sake of completeness, and "
          "have not been fully verified.");
      debugOutput(message, message, COLOR_YELLOW);

      if (_is_eigen && _exec_type == ExecutionType::Transient)
        mooseError(
            "Eigenvalue kernels are not compatible with transient executioners. Did you intend "
            "to run a steady-state simulation?");
    }

    if (_is_eigen && getParam<bool>("debug_disable_fission") && _particle == Particletype::Neutron)
      mooseError("Fission cannot be disabled in eigenvalue simulations.");
  }

  // Initial the diffusion approximation separately.
  if (_transport_scheme == TransportScheme::DiffusionApprox)
  {
    if (!_var_init && _problem)
    {
      debugOutput("  - Scheme: Diffusion-CGFEM", "  - Scheme: Diffusion-CGFEM");

      // Set the enum so quadrature sets can be determined appropriately.
      for (unsigned int g = 0; g < _num_groups; ++g)
      {
        // Set up variable names for the group flux moments.
        _group_flux_moments.emplace(g, std::vector<VariableName>());
        _group_flux_moments[g].emplace_back(_flux_moment_name + "_" + Moose::stringify(g + 1u) +
                                            "_" + Moose::stringify(0u) + "_" +
                                            Moose::stringify(0u));
      }

      _var_init = true;
      debugOutput("  - Initializing flux groups...", "  - Initializing flux groups...");
    }
  }
  else if (_transport_scheme == TransportScheme::FluxMomentTransfer)
  {
    if (!_var_init && _problem)
    {
      debugOutput("  - Scheme: Flux Moment Transfer", "  - Scheme: Flux Moment Transfer");

      // Set the enum so quadrature sets can be determined appropriately.
      for (unsigned int g = 0; g < _num_groups; ++g)
      {
        // Set up variable names for the group flux moments.
        _group_flux_moments.emplace(g, std::vector<VariableName>());
        _group_flux_moments[g].emplace_back(_flux_moment_name + "_" + Moose::stringify(g + 1u) +
                                            "_" + Moose::stringify(0u) + "_" +
                                            Moose::stringify(0u));
      }

      _var_init = true;
      debugOutput("  - Initializing flux groups...", "  - Initializing flux groups...");
    }
  }
  else
  {
    if (!_var_init && _problem)
    {
      debugOutput("Transport System Initialization: ", "Transport System Initialization: ");
      debugOutput("  - Scheme: SAAF-CGFEM", "  - Scheme: SAAF-CGFEM");
      debugOutput("  - Building the angular quadrature set...",
                  "  - Building the angular quadrature set...");

      // Set the enum so quadrature sets can be determined appropriately.
      _num_group_moments = 0u;
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
              _group_flux_moments[g].emplace_back(_flux_moment_name + "_" +
                                                  Moose::stringify(g + 1u) + "_" +
                                                  Moose::stringify(l) + "_" + Moose::stringify(0));
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
                _group_flux_moments[g].emplace_back(
                    _flux_moment_name + "_" + Moose::stringify(g + 1u) + "_" + Moose::stringify(l) +
                    "_" + Moose::stringify(m));
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
                _group_flux_moments[g].emplace_back(
                    _flux_moment_name + "_" + Moose::stringify(g + 1u) + "_" + Moose::stringify(l) +
                    "_" + Moose::stringify(m));
              }
            }
          }
          break;

        default:
          mooseError("Unknown mesh dimensionality.");
          break;
      }

      switch (_p_type)
      {
        case ProblemType::Cartesian1D:
          _num_flux_ordinates = 2u * _n_l;
          break;

        case ProblemType::Cartesian2D:
          _num_flux_ordinates = 4u * _n_l * _n_c;
          break;

        case ProblemType::Cartesian3D:
          _num_flux_ordinates = 8u * _n_l * _n_c;
          break;

        default:
          mooseError("Unknown mesh dimensionality.");
          break;
      }

      std::string set_info(
          std::string("    - Angular quadrature set information:\n") +
          std::string("      - Polar points per quadrant: ") + Moose::stringify(_n_l) +
          std::string("\n      - Azimuthal points per quadrant: ") + Moose::stringify(_n_c) +
          std::string("\n      - Total number of flux ordinates: ") +
          Moose::stringify(_num_flux_ordinates));
      debugOutput(set_info);

      debugOutput("  - Determining variable names...", "  - Determining variable names...");
      for (unsigned int g = 0; g < _num_groups; ++g)
      {
        // Set up variable names for the group angular fluxes.
        _group_angular_fluxes.emplace(g, std::vector<VariableName>());
        _group_angular_fluxes[g].reserve(_num_flux_ordinates);
        for (unsigned int n = 1; n <= _num_flux_ordinates; ++n)
        {
          _group_angular_fluxes[g].emplace_back(
              _angular_flux_name + "_" + Moose::stringify(g + 1u) + "_" + Moose::stringify(n));
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

          addSNBCs(var_name, g, n);
        }

        // Add initial conditions.
        if (_current_task == "add_ic")
        {
          if (g == 0u && n == 0u)
            debugOutput("    - Adding ICs...");

          addSNICs(var_name, g);
        }
      }

      // Loop over all moments and set up auxvariables and auxkernels.
      unsigned int moment_index = 0u;
      switch (_p_type)
      {
        case ProblemType::Cartesian1D:
          for (unsigned int l = 0; l <= _max_eval_anisotropy; ++l)
          {
            const auto & var_name = _group_flux_moments[g][moment_index];

            // Add auxvariables.
            if (_current_task == "add_aux_variable")
            {
              if (g == 0u && l == 0u)
                debugOutput("    - Adding auxvariables...");

              addAuxVariables(var_name);
            }

            // Add auxkernels.
            if (_current_task == "add_aux_kernel")
            {
              if (g == 0u && l == 0u)
                debugOutput("    - Adding auxkernels...");

              addAuxKernels(var_name, g, l, 0u);
            }

            moment_index++;
          }
          moment_index = 0u;
          break;

        case ProblemType::Cartesian2D:
          for (unsigned int l = 0u; l <= _max_eval_anisotropy; ++l)
          {
            for (int m = 0; m <= static_cast<int>(l); ++m)
            {
              const auto & var_name = _group_flux_moments[g][moment_index];

              // Add auxvariables.
              if (_current_task == "add_aux_variable")
              {
                if (g == 0u && l == 0u && m == 0u)
                  debugOutput("    - Adding auxvariables...");

                addAuxVariables(var_name);
              }

              // Add auxkernels.
              if (_current_task == "add_aux_kernel")
              {
                if (g == 0u && l == 0u && m == 0u)
                  debugOutput("    - Adding auxkernels...");

                addAuxKernels(var_name, g, l, m);
              }

              moment_index++;
            }
          }
          moment_index = 0u;
          break;

        case ProblemType::Cartesian3D:
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
          moment_index = 0u;
          break;

        default:
          mooseError("Unknown mesh dimensionality.");
          break;
      }
    }
  }

  // Inform the system that it needs to restart.
  if (_current_task == "check_copy_nodal_vars" && getParam<bool>("init_from_file"))
    _app.setExodusFileRestart(true);

  // Add all required variables to a list for reinitialization from a file.
  if (_current_task == "copy_nodal_vars" && getParam<bool>("init_from_file"))
  {
    auto & system = _problem->getNonlinearSystemBase();
    auto & aux_system = _problem->getAuxiliarySystem();

    switch (_transport_scheme)
    {
      case TransportScheme::SAAFCFEM:
        for (unsigned int g = 0u; g < _num_groups; ++g)
        {
          for (unsigned int n = 0u; n < _num_flux_ordinates; ++n)
            system.addVariableToCopy(
                _group_angular_fluxes[g][n], _group_angular_fluxes[g][n], "LATEST");

          for (unsigned int m = 0u; m < _num_group_moments; ++m)
            aux_system.addVariableToCopy(
                _group_flux_moments[g][m], _group_flux_moments[g][m], "LATEST");
        }
        break;
      case TransportScheme::DiffusionApprox:
        for (unsigned int g = 0u; g < _num_groups; ++g)
          system.addVariableToCopy(_group_flux_moments[g][0], _group_flux_moments[g][0], "LATEST");

        break;
      case TransportScheme::FluxMomentTransfer:
        for (unsigned int g = 0u; g < _num_groups; ++g)
          aux_system.addVariableToCopy(
              _group_flux_moments[g][0], _group_flux_moments[g][0], "LATEST");

        break;
    }
  }

  if (_current_task == "add_postprocessor" && getParam<bool>("is_conservative_transfer_src"))
  {
    debugOutput("    - Add post-processors...");

    for (unsigned int g = 0u; g < _num_groups; ++g)
      addSourceConservativePP(_group_flux_moments[g][0]);
  }

  // This should really be it's own thing on the task graph instead of relying on 'add_variable'
  // being after 'common_output'. TODO: change this.
  if (_current_task == "add_variable")
  {
    debugOutput("    - Modifying outputs...");

    if (!getParam<bool>("output_angular_fluxes") && _transport_scheme == TransportScheme::SAAFCFEM)
      modifyOutputs();
  }

  // Add SN user objects, namely the quadrature set.
  if (_current_task == "add_user_object" && _transport_scheme == TransportScheme::SAAFCFEM)
  {
    debugOutput("    - Add user objects...");

    addSNUserObjects();
  }
}
//------------------------------------------------------------------------------

//------------------------------------------------------------------------------
// Initialize the CGFEM-SAAF scheme.
//------------------------------------------------------------------------------
void
TransportAction::actSAAFCFEM()
{
  // Loop over all groups.
  for (unsigned int g = 0; g < _num_groups; ++g)
  {
    // Loop over all the flux ordinates.
    for (unsigned int n = 0; n < _num_flux_ordinates; ++n)
    {
      const auto & var_name = _group_angular_fluxes[g][n];

      if (_current_task == "add_dirac_kernel")
      {
        if (g == 0u && n == 0u)
          debugOutput("    - Adding Dirac kernels...");

        addSAAFDiracKernels(var_name, g, n);
      }

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
}
//------------------------------------------------------------------------------

//------------------------------------------------------------------------------
// Initialize the diffusion approximation scheme.
//------------------------------------------------------------------------------
void
TransportAction::actDiffusion()
{
  // Loop over all groups.
  for (unsigned int g = 0; g < _num_groups; ++g)
  {
    const auto & var_name = _group_flux_moments[g][0];

    // Add a non-linear variable.
    if (_current_task == "add_variable")
    {
      if (g == 0u)
        debugOutput("    - Adding variables...");

      addVariable(var_name);
    }

    // Add boundary conditions.
    if (_current_task == "add_bc")
    {
      if (g == 0u)
        debugOutput("    - Adding BCs...");

      addDiffusionBCs(var_name);
    }

    // Add initial conditions.
    if (_current_task == "add_ic")
    {
      if (g == 0u)
        debugOutput("    - Adding ICs...");

      addDiffusionICs(var_name, g);
    }

    // Add kernels.
    if (_current_task == "add_kernel")
    {
      if (g == 0u)
        debugOutput("    - Adding kernels...");

      addDiffusionKernels(var_name, g);
    }

    // Add Dirac kernels.
    if (_current_task == "add_dirac_kernel")
    {
      if (g == 0u)
        debugOutput("    - Adding Dirac kernels...");

      addDiffusionDiracKernels(var_name, g);
    }
  }
}
//------------------------------------------------------------------------------

//------------------------------------------------------------------------------
// Initialize the transfer scheme. Used to pull in flux moments from a sub-app in a
// multi-app simulation.
//------------------------------------------------------------------------------
void
TransportAction::actTransfer()
{
  std::string source_var_name = "";
  // Loop over all groups.
  for (unsigned int g = 0; g < _num_groups; ++g)
  {
    const auto & var_name = _group_flux_moments[g][0];

    source_var_name = _source_flux_moment_names + "_" + Moose::stringify(g + 1u) + "_0_0";

    // Add the resulting transferred variable.
    if (_current_task == "add_aux_variable")
    {
      if (g == 0u)
        debugOutput("    - Adding auxvariables...");

      addAuxVariables(var_name);
    }

    // Add a transfer.
    if (_current_task == "add_transfer" && !getParam<bool>("init_from_file"))
    {
      if (g == 0u)
        debugOutput("    - Adding transfers...");

      addTransfers(var_name, source_var_name);
    }

    if (_current_task == "add_postprocessor" && !getParam<bool>("init_from_file"))
      addDestinationConservativePP(var_name);
  }
}
//------------------------------------------------------------------------------

void
TransportAction::actUncollided()
{
  std::string unc_source_var_name = "";
  std::string unc_var_name = "";
  // Loop over all groups.
  for (unsigned int g = 0; g < _num_groups; ++g)
  {
    unc_var_name = _flux_moment_name + "_" + Moose::stringify(g + 1u) + "_0_0_uncollided";
    unc_source_var_name =
        _uncollided_source_flux_moment_names + "_" + Moose::stringify(g + 1u) + "_0_0";

    // Add the resulting transferred variable.
    if (_current_task == "add_aux_variable")
    {
      if (g == 0u)
        debugOutput("    - Adding auxvariables...");

      addAuxVariables(unc_var_name);
    }

    // Add a transfer.
    if (_current_task == "add_transfer" && !getParam<bool>("init_from_file"))
    {
      if (g == 0u)
        debugOutput("    - Adding transfers...");

      addUncTransfers(unc_var_name, unc_source_var_name);
    }
  }
}

//------------------------------------------------------------------------------
// Functions to add transfers for decoupled uncollided flux calculations.
//------------------------------------------------------------------------------
void
TransportAction::addUncTransfers(const std::string & to_var_name,
                                 const std::string & source_var_name)
{
  auto params = _factory.getValidParams("MultiAppGeneralFieldShapeEvaluationTransfer");
  params.set<std::vector<VariableName>>("source_variable").emplace_back(source_var_name);
  params.set<std::vector<AuxVariableName>>("variable").emplace_back(to_var_name);

  params.set<MultiAppName>("from_multi_app") = _uncollided_from_multi_app_name;

  auto & from_blocks = params.set<std::vector<SubdomainName>>("from_blocks");
  for (const auto & block : _uncollided_from_subdomain_ids)
    from_blocks.emplace_back(block);

  auto & to_block = params.set<std::vector<SubdomainName>>("to_blocks");
  for (const auto & block : getParam<std::vector<SubdomainName>>("block"))
    to_block.emplace_back(block);

  /*
  if (getParam<bool>("use_conservative_transfers"))
  {
    params.set<std::vector<PostprocessorName>>("from_postprocessors_to_be_preserved")
        .emplace_back("ElementIntegralVariablePostprocessor_" + source_var_name + "_src");
    params.set<std::vector<PostprocessorName>>("to_postprocessors_to_be_preserved")
        .emplace_back("ElementIntegralVariablePostprocessor_" + to_var_name + "_dst");
  }
  */

  _problem->addTransfer("MultiAppGeneralFieldShapeEvaluationTransfer",
                        "MultiAppGeneralFieldShapeEvaluationTransfer_" + to_var_name + "_from_" +
                            source_var_name,
                        params);
  debugOutput(
      "      - Adding Transfer MultiAppGeneralFieldShapeEvaluationTransfer for the variable " +
      to_var_name + " from " + source_var_name + ".");
}

//------------------------------------------------------------------------------
// Functions to add transfers for decoupled multi-app calculations.
//------------------------------------------------------------------------------
void
TransportAction::addTransfers(const std::string & to_var_name, const std::string & source_var_name)
{
  if (!getParam<bool>("use_copy"))
  {
    auto params = _factory.getValidParams("MultiAppGeneralFieldShapeEvaluationTransfer");
    params.set<std::vector<VariableName>>("source_variable").emplace_back(source_var_name);
    params.set<std::vector<AuxVariableName>>("variable").emplace_back(to_var_name);

    if (_from_multi_app_name != "")
      params.set<MultiAppName>("from_multi_app") = _from_multi_app_name;

    auto & from_blocks = params.set<std::vector<SubdomainName>>("from_blocks");
    for (const auto & block : _from_subdomain_ids)
      from_blocks.emplace_back(block);

    auto & to_block = params.set<std::vector<SubdomainName>>("to_blocks");
    for (const auto & block : getParam<std::vector<SubdomainName>>("block"))
      to_block.emplace_back(block);

    if (getParam<bool>("use_conservative_transfers"))
    {
      params.set<std::vector<PostprocessorName>>("from_postprocessors_to_be_preserved")
          .emplace_back("ElementIntegralVariablePostprocessor_" + source_var_name + "_src");
      params.set<std::vector<PostprocessorName>>("to_postprocessors_to_be_preserved")
          .emplace_back("ElementIntegralVariablePostprocessor_" + to_var_name + "_dst");
    }

    _problem->addTransfer("MultiAppGeneralFieldShapeEvaluationTransfer",
                          "MultiAppGeneralFieldShapeEvaluationTransfer_" + to_var_name + "_from_" +
                              source_var_name,
                          params);
    debugOutput(
        "      - Adding Transfer MultiAppGeneralFieldShapeEvaluationTransfer for the variable " +
        to_var_name + " from " + source_var_name + ".");
  }
  else
  {
    auto params = _factory.getValidParams("MultiAppCopyTransfer");
    params.set<std::vector<VariableName>>("source_variable").emplace_back(source_var_name);
    params.set<std::vector<AuxVariableName>>("variable").emplace_back(to_var_name);

    if (_from_multi_app_name != "")
      params.set<MultiAppName>("from_multi_app") = _from_multi_app_name;

    _problem->addTransfer("MultiAppCopyTransfer",
                          "MultiAppCopyTransfer_" + to_var_name + "_from_" + source_var_name,
                          params);
    debugOutput("      - Adding Transfer MultiAppCopyTransfer for the variable " + to_var_name +
                " from " + source_var_name + ".");
  }
}

void
TransportAction::addDestinationConservativePP(const std::string & to_var_name)
{
  auto params = _factory.getValidParams("ElementIntegralVariablePostprocessor");
  params.set<std::vector<VariableName>>("variable").emplace_back(to_var_name);
  params.set<ExecFlagEnum>("execute_on") = EXEC_TRANSFER;
  params.set<std::vector<OutputName>>("outputs").emplace_back("none");

  if (isParamValid("block"))
    params.set<std::vector<SubdomainName>>("block") = getParam<std::vector<SubdomainName>>("block");

  _problem->addPostprocessor("ElementIntegralVariablePostprocessor",
                             "ElementIntegralVariablePostprocessor_" + to_var_name + "_dst",
                             params);
  debugOutput("      - Adding Transfer ElementIntegralVariablePostprocessor for the variable " +
              to_var_name + ".");
}

void
TransportAction::addSourceConservativePP(const std::string & source_var_name)
{
  auto params = _factory.getValidParams("ElementIntegralVariablePostprocessor");
  params.set<std::vector<VariableName>>("variable").emplace_back(source_var_name);
  params.set<ExecFlagEnum>("execute_on") = EXEC_TIMESTEP_END;
  params.set<std::vector<OutputName>>("outputs").emplace_back("none");

  if (isParamValid("block"))
    params.set<std::vector<SubdomainName>>("block") = getParam<std::vector<SubdomainName>>("block");

  _problem->addPostprocessor("ElementIntegralVariablePostprocessor",
                             "ElementIntegralVariablePostprocessor_" + source_var_name + "_src",
                             params);
  debugOutput("      - Adding Transfer ElementIntegralVariablePostprocessor for the variable " +
              source_var_name + ".");
}
//------------------------------------------------------------------------------

//------------------------------------------------------------------------------
// Functions to add MOOSE objects common to all schemes.
//------------------------------------------------------------------------------
void
TransportAction::modifyOutputs()
{
  // Fetch all AddOutputAction's from the action warehouse.
  const auto & output_actions = _app.actionWarehouse().getActionListByName("add_output");
  for (const auto & act : output_actions)
  {
    // Extract the Output action.
    AddOutputAction * action = dynamic_cast<AddOutputAction *>(act);
    if (!action)
      continue;

    InputParameters & output_params = action->getObjectParams();
    if (output_params.have_parameter<std::vector<VariableName>>("hide"))
    {
      for (const auto & [g, flux_ordinate] : _group_angular_fluxes)
      {
        for (unsigned int n = 0u; n < flux_ordinate.size(); ++n)
          output_params.set<std::vector<VariableName>>("hide").emplace_back(flux_ordinate[n]);
      }
    }
  }
}

void
TransportAction::addSNUserObjects()
{
  // Add AQProvider.
  {
    auto params = _factory.getValidParams("AQProvider");

    // Assign dimensionality.
    MooseEnum dimensionality("1D_cartesian 2D_cartesian 3D_cartesian");
    dimensionality.assign(static_cast<int>(_p_type));
    params.set<MooseEnum>("dimensionality") = dimensionality;

    // Assign major axis.
    params.set<MooseEnum>("major_axis") = getParam<MooseEnum>("major_axis");

    // Assign quadrature rules.
    params.set<unsigned int>("n_l") = 2u * _n_l;
    params.set<unsigned int>("n_c") = 2u * _n_c;

    _problem->addUserObject("AQProvider", "AQProvider_" + name(), params);
    debugOutput("      - Adding UserObject AQProvider_" + name() + ".");
  } // AQProvider
}

void
TransportAction::addSNBCs(const std::string & var_name, unsigned int g, unsigned int n)
{
  // Add SNVacuumBC.
  if (_vacuum_side_sets.size() > 0u)
  {
    auto params = _factory.getValidParams("SNVacuumBC");
    params.set<NonlinearVariableName>("variable") = var_name;
    // Ordinate index is required to fetch the particle direction.
    params.set<unsigned int>("ordinate_index") = n;

    // Apply the parameters for the quadrature rule.
    applyQuadratureParameters(params);

    params.set<std::vector<BoundaryName>>("boundary") = _vacuum_side_sets;

    _problem->addBoundaryCondition("SNVacuumBC", "SNVacuumBC_" + var_name, params);
    debugOutput("      - Adding BC SNVacuumBC for the variable " + var_name + ".");
  } // SNVacuumBC

  // Add SNSourceBC.
  if (_source_side_sets.size() > 0u && !_using_uncollided)
  {
    if (_source_side_sets.size() != _boundary_source_moments.size() &&
        _source_side_sets.size() != _boundary_source_anisotropy.size())
    {
      mooseError("There is a mismatch between the number of source boundary conditions and the "
                 "number of provided moments / anisotropy.");
    }

    for (unsigned int i = 0u; i < _source_side_sets.size(); ++i)
    {
      auto params = _factory.getValidParams("SNSourceBC");
      params.set<NonlinearVariableName>("variable") = var_name;
      params.set<unsigned int>("num_groups") = _num_groups;
      params.set<unsigned int>("group_index") = g;
      // Ordinate index is required to fetch the particle direction.
      params.set<unsigned int>("ordinate_index") = n;

      // The group source and it's degree of anisotropy.
      for (unsigned int m = 0u; m < _boundary_source_moments[i].size(); ++m)
        params.set<std::vector<Real>>("group_source").emplace_back(_boundary_source_moments[i][m]);
      params.set<unsigned int>("source_anisotropy") = _boundary_source_anisotropy[i];

      // Apply the parameters for the quadrature rule.
      applyQuadratureParameters(params);

      params.set<std::vector<BoundaryName>>("boundary").emplace_back(_source_side_sets[i]);

      _problem->addBoundaryCondition("SNSourceBC",
                                     "SNSourceBC_" + var_name + "_" +
                                         Moose::stringify(_source_side_sets[i]),
                                     params);
      debugOutput("Adding BC SNSourceBC for the variable " + var_name + ".");
    }
  } // SNSourceBC

  // Add SNNormalCurrentBC.
  if (_current_side_sets.size() > 0u && !_using_uncollided)
  {
    if (_current_side_sets.size() != _boundary_currents.size() &&
        _current_side_sets.size() != _boundary_current_anisotropy.size())
    {
      mooseError("There is a mismatch between the number of current boundary conditions and the "
                 "number of provided currents / the approximation anisotropy.");
    }

    for (unsigned int i = 0u; i < _current_side_sets.size(); ++i)
    {
      auto params = _factory.getValidParams("SNNormalCurrentBC");
      params.set<NonlinearVariableName>("variable") = var_name;
      params.set<unsigned int>("num_groups") = _num_groups;
      params.set<unsigned int>("group_index") = g;
      // Ordinate index is required to fetch the particle direction.
      params.set<unsigned int>("ordinate_index") = n;

      // The group currents and the degree of anisotropy to apply to the approximation.
      params.set<Real>("group_current") = _boundary_currents[i][g];
      params.set<unsigned int>("current_anisotropy") = _boundary_current_anisotropy[i];

      // Apply the parameters for the quadrature rule.
      applyQuadratureParameters(params);

      params.set<std::vector<BoundaryName>>("boundary").emplace_back(_current_side_sets[i]);

      _problem->addBoundaryCondition("SNNormalCurrentBC",
                                     "SNNormalCurrentBC_" + var_name + "_" +
                                         Moose::stringify(_current_side_sets[i]),
                                     params);
      debugOutput("Adding BC SNNormalCurrentBC for the variable " + var_name + ".");
    }
  } // SNNormalCurrentBC

  // Add ADSNReflectiveBC.
  if (_reflective_side_sets.size() > 0u)
  {
    auto params = _factory.getValidParams("ADSNReflectiveBC");
    params.set<NonlinearVariableName>("variable") = var_name;
    // Ordinate index is required to fetch the particle direction.
    params.set<unsigned int>("ordinate_index") = n;

    // The flux ordinates for this group.
    params.set<std::vector<VariableName>>("psi_ref") = _group_angular_fluxes[g];

    // Apply the parameters for the quadrature rule.
    applyQuadratureParameters(params);

    params.set<std::vector<BoundaryName>>("boundary") = _reflective_side_sets;

    _problem->addBoundaryCondition("ADSNReflectiveBC", "ADSNReflectiveBC_" + var_name, params);
    debugOutput("Adding BC ADSNReflectiveBC for the variable " + var_name + ".");
  } // ADSNReflectiveBC
}

// TODO: Initial conditions.
void
TransportAction::addSNICs(const std::string & var_name, unsigned int g)
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

        const auto & const_ic = getParam<std::vector<Real>>("constant_ic");
        if (const_ic.size() == 1)
          params.set<Real>("value") = const_ic[0];
        else if (const_ic.size() == _num_groups)
          params.set<Real>("value") = const_ic[g];
        else
          mooseError("Size of 'constant_ic' does not match the declared number "
                     "of particle groups.");

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
TransportAction::addAuxVariables(const std::string & var_name)
{
  if (getParam<bool>("transfer_to_fv") && _transport_scheme == TransportScheme::FluxMomentTransfer)
  {
    auto fe_type = AddVariableAction::feType(_pars);
    auto type = AddVariableAction::variableType(fe_type, true, false);
    auto params = _factory.getValidParams(type);
    params.set<MooseEnum>("order") = fe_type.order.get_order();
    params.set<MooseEnum>("family") = Moose::stringify(fe_type.family);

    if (isParamValid("block"))
    {
      params.set<std::vector<SubdomainName>>("block") =
          getParam<std::vector<SubdomainName>>("block");
    }

    _problem->addAuxVariable(type, var_name, params);
    debugOutput("      - Adding auxvariable " + var_name + ".");
  }
  else
  {
    auto fe_type = AddVariableAction::feType(_pars);
    auto type = AddVariableAction::variableType(fe_type, false, false);
    auto params = _factory.getValidParams(type);
    params.set<MooseEnum>("order") = fe_type.order.get_order();
    params.set<MooseEnum>("family") = Moose::stringify(fe_type.family);

    if (isParamValid("block"))
    {
      params.set<std::vector<SubdomainName>>("block") =
          getParam<std::vector<SubdomainName>>("block");
    }

    _problem->addAuxVariable(type, var_name, params);
    debugOutput("      - Adding auxvariable " + var_name + ".");
  }
}

void
TransportAction::addAuxKernels(const std::string & var_name, unsigned int g, unsigned int l, int m)
{
  // Add ParticleFluxMoment.
  {
    InputParameters params = _factory.getValidParams("ParticleFluxMoment");
    params.set<AuxVariableName>("variable") = var_name;
    // Flux moment degree and order.
    params.set<unsigned int>("degree") = l;
    params.set<int>("order") = m;
    params.set<unsigned int>("group_index") = g;
    params.set<unsigned int>("num_groups") = _num_groups;

    params.set<ExecFlagEnum>("execute_on") = {EXEC_INITIAL, EXEC_TIMESTEP_BEGIN, EXEC_LINEAR};

    // Apply the parameters for the quadrature rule.
    applyQuadratureParameters(params);

    // The flux ordinates for this group.
    params.set<std::vector<VariableName>>("group_flux_ordinates") = _group_angular_fluxes[g];

    if (_using_uncollided && l == 0u && m == 0)
      params.set<std::vector<VariableName>>("uncollided_flux_moment")
          .emplace_back(_flux_moment_name + "_" + Moose::stringify(g + 1u) + "_0_0_uncollided");

    if (isParamValid("block"))
    {
      params.set<std::vector<SubdomainName>>("block") =
          getParam<std::vector<SubdomainName>>("block");
    }

    _problem->addAuxKernel("ParticleFluxMoment", "ParticleFluxMoment_" + var_name, params);
    debugOutput("      - Adding auxkernel ParticleFluxMoment for the variable " + var_name + ".");
  } // ParticleFluxMoment

  // Add a ParticleFluxMoment which scales the moments at the end of the solve..
  if (getParam<bool>("scale_sources"))
  {
    InputParameters params = _factory.getValidParams("ParticleFluxMoment");
    params.set<AuxVariableName>("variable") = var_name;
    // Flux moment degree and order.
    params.set<unsigned int>("degree") = l;
    params.set<int>("order") = m;
    params.set<unsigned int>("group_index") = g;
    params.set<unsigned int>("num_groups") = _num_groups;

    if (!_using_uncollided)
      params.set<Real>("scale_factor") = _source_scale_factor;

    params.set<ExecFlagEnum>("execute_on") = {EXEC_TIMESTEP_END};

    // Apply the parameters for the quadrature rule.
    applyQuadratureParameters(params);

    // The flux ordinates for this group.
    params.set<std::vector<VariableName>>("group_flux_ordinates") = _group_angular_fluxes[g];

    if (_using_uncollided && l == 0u && m == 0)
      params.set<std::vector<VariableName>>("uncollided_flux_moment")
          .emplace_back(_flux_moment_name + "_" + Moose::stringify(g + 1u) + "_0_0_uncollided");

    if (isParamValid("block"))
    {
      params.set<std::vector<SubdomainName>>("block") =
          getParam<std::vector<SubdomainName>>("block");
    }

    _problem->addAuxKernel("ParticleFluxMoment", "Scaling_ParticleFluxMoment_" + var_name, params);
    debugOutput("      - Adding auxkernel ParticleFluxMoment to scale the variable " + var_name +
                ".");
  } // ParticleFluxMoment
}
//------------------------------------------------------------------------------

//------------------------------------------------------------------------------
// Functions to add MOOSE objects for the CGFEM-SAAF scheme.
//------------------------------------------------------------------------------
void
TransportAction::addSAAFKernels(const std::string & var_name, unsigned int g, unsigned int n)
{
  // Add SAAFTimeDerivative.
  if (_exec_type == ExecutionType::Transient)
  {
    auto params = _factory.getValidParams("SAAFTimeDerivative");
    params.set<NonlinearVariableName>("variable") = var_name;
    // Set the name of the TransportAction so it can fetch the appropriate material properties.
    params.set<std::string>("transport_system") = name();
    // Group index is required to fetch the group particle velocity.
    params.set<unsigned int>("group_index") = g;
    // Ordinate index is required to fetch the direction of travel for
    // stabilization.
    params.set<unsigned int>("ordinate_index") = n;

    // Apply the parameters for the quadrature rule.
    applyQuadratureParameters(params);

    if (isParamValid("block"))
    {
      params.set<std::vector<SubdomainName>>("block") =
          getParam<std::vector<SubdomainName>>("block");
    }

    _problem->addKernel("SAAFTimeDerivative", "SAAFTimeDerivative_" + var_name, params);
    debugOutput("      - Adding kernel SAAFTimeDerivative for the variable " + var_name + ".");
  } // SAAFTimeDerivative

  // Add SAAFStreaming.
  {
    auto params = _factory.getValidParams("SAAFStreaming");
    params.set<NonlinearVariableName>("variable") = var_name;
    // Set the name of the TransportAction so it can fetch the appropriate material properties.
    params.set<std::string>("transport_system") = name();
    // Group index is required to fetch the group particle removal cross-section
    // for stabilization.
    params.set<unsigned int>("group_index") = g;
    // Ordinate index is required to fetch the particle direction.
    params.set<unsigned int>("ordinate_index") = n;

    // Apply the parameters for the quadrature rule.
    applyQuadratureParameters(params);

    if (isParamValid("block"))
    {
      params.set<std::vector<SubdomainName>>("block") =
          getParam<std::vector<SubdomainName>>("block");
    }

    _problem->addKernel("SAAFStreaming", "SAAFStreaming_" + var_name, params);
    debugOutput("      - Adding kernel SAAFStreaming for the variable " + var_name + ".");
  } // SAAFStreaming

  // Add SNRemoval.
  {
    auto params = _factory.getValidParams("SNRemoval");
    params.set<NonlinearVariableName>("variable") = var_name;
    // Set the name of the TransportAction so it can fetch the appropriate material properties.
    params.set<std::string>("transport_system") = name();
    // Group index is required to fetch the group particle removal
    // cross-section.
    params.set<unsigned int>("group_index") = g;

    if (isParamValid("block"))
    {
      params.set<std::vector<SubdomainName>>("block") =
          getParam<std::vector<SubdomainName>>("block");
    }

    _problem->addKernel("SNRemoval", "SNRemoval_" + var_name, params);
    debugOutput("      - Adding kernel SNRemoval for the variable " + var_name + ".");
  } // SNRemoval

  // Add SAAFVolumeSource
  if (_volumetric_source_blocks.size() > 0u && !_is_eigen && !_using_uncollided)
  {
    for (unsigned int i = 0u; i < _volumetric_source_blocks.size(); ++i)
    {
      auto params = _factory.getValidParams("SAAFVolumeSource");
      params.set<NonlinearVariableName>("variable") = var_name;
      // Set the name of the TransportAction so it can fetch the appropriate material properties.
      params.set<std::string>("transport_system") = name();
      // Group index is required to fetch the group particle removal cross-section
      // for stabilization.
      params.set<unsigned int>("group_index") = g;
      // Number of groups
      params.set<unsigned int>("num_groups") = _num_groups;
      // Ordinate index is required to fetch the particle direction.
      params.set<unsigned int>("ordinate_index") = n;

      params.set<std::vector<Real>>("group_source") = _volumetric_source_moments[i];
      params.set<unsigned int>("source_anisotropy") = _volumetric_source_anisotropy[i];

      // Apply the parameters for the quadrature rule.
      applyQuadratureParameters(params);

      params.set<std::vector<SubdomainName>>("block").emplace_back(_volumetric_source_blocks[i]);

      _problem->addKernel("SAAFVolumeSource",
                          "SAAFVolumeSource_" + var_name + "_" + _volumetric_source_blocks[i],
                          params);
      debugOutput("      - Adding kernel SAAFVolumeSource for the variable " + var_name + ".");
    }
  } // SAAFVolumeSource

  // Add SAAFMaterialSource.
  if (!_is_eigen && !_using_uncollided)
  {
    auto params = _factory.getValidParams("SAAFMaterialSource");
    params.set<NonlinearVariableName>("variable") = var_name;
    // Set the name of the TransportAction so it can fetch the appropriate material properties.
    params.set<std::string>("transport_system") = name();
    // Group index and the number of groups are required to fetch the
    // source moments for the proper spectral energy group.
    params.set<unsigned int>("group_index") = g;
    params.set<unsigned int>("num_groups") = _num_groups;
    // Ordinate index is required to fetch the particle direction.
    params.set<unsigned int>("ordinate_index") = n;

    // Apply the parameters for the quadrature rule.
    applyQuadratureParameters(params);

    if (isParamValid("block"))
    {
      params.set<std::vector<SubdomainName>>("block") =
          getParam<std::vector<SubdomainName>>("block");
    }

    _problem->addKernel("SAAFMaterialSource", "SAAFMaterialSource_" + var_name, params);
    debugOutput("      - Adding kernel SAAFMaterialSource for the variable " + var_name + ".");
  } // SAAFMaterialSource

  // Only add fission kernels if debug doesn't disable them AND this transport system represents a
  // neutron field.
  if (!getParam<bool>("debug_disable_fission") && _particle == Particletype::Neutron)
  {
    // Add SAAFMomentFission.
    {
      auto params = _factory.getValidParams("SAAFMomentFission");
      params.set<NonlinearVariableName>("variable") = var_name;
      // Set the name of the TransportAction so it can fetch the appropriate material
      // properties.
      params.set<std::string>("transport_system") = name();
      // Group index and the number of groups are required to fetch the
      // scattering cross-section moments.
      params.set<unsigned int>("group_index") = g;
      params.set<unsigned int>("num_groups") = _num_groups;
      // Ordinate index is required to fetch the particle direction.
      params.set<unsigned int>("ordinate_index") = n;

      // Apply the parameters for the quadrature rule.
      applyQuadratureParameters(params);

      // Copy all of the scalar flux names into the variable parameter.
      auto & scalar_flux_names = params.set<std::vector<VariableName>>("group_scalar_fluxes");
      for (unsigned int g_prime = 0; g_prime < _num_groups; ++g_prime)
        scalar_flux_names.emplace_back(_group_flux_moments[g_prime][0u]);

      if (isParamValid("block"))
      {
        params.set<std::vector<SubdomainName>>("block") =
            getParam<std::vector<SubdomainName>>("block");
      }

      // For eigenvalues.
      if (_is_eigen)
        params.set<std::vector<TagName>>("extra_vector_tags").emplace_back("eigen");

      _problem->addKernel("SAAFMomentFission", "SAAFMomentFission_" + var_name, params);
      debugOutput("      - Adding kernel SAAFMomentFission for the variable " + var_name + ".");
    } // SAAFMomentFission
  }

  // Only add scattering kernels if debug doesn't disable them.
  if (!getParam<bool>("debug_disable_scattering"))
  {
    // Debug option to disable the source iteration solver.
    if (!getParam<bool>("debug_disable_source_iteration"))
    {
      // Compute the scattering evaluation with source iteration.
      mooseError("Scattering iteration is not supported and is a work in "
                 "progress.");
    }
    else
    {
      if (getParam<bool>("use_scattering_jacobians"))
      {
        // Computes the scattering evaluation without source iteration using a hand-coded Jacobian.
        // Add SAAFScattering.
        {
          auto params = _factory.getValidParams("SAAFScattering");
          params.set<NonlinearVariableName>("variable") = var_name;
          // Set the name of the TransportAction so it can fetch the appropriate material
          // properties.
          params.set<std::string>("transport_system") = name();
          // Group index and the number of groups are required to fetch the
          // scattering cross-section moments.
          params.set<unsigned int>("group_index") = g;
          params.set<unsigned int>("num_groups") = _num_groups;
          // Ordinate index is required to fetch the particle direction.
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
            params.set<std::vector<SubdomainName>>("block") =
                getParam<std::vector<SubdomainName>>("block");
          }

          _problem->addKernel("SAAFScattering", "SAAFScattering_" + var_name, params);
          debugOutput("      - Adding kernel SAAFScattering for the variable " + var_name + ".");
        } // SAAFScattering
      }
      else
      {
        // Computes the scattering evaluation without source iteration without a Jacobian.
        // Add SAAFMomentScattering.
        {
          auto params = _factory.getValidParams("SAAFMomentScattering");
          params.set<NonlinearVariableName>("variable") = var_name;
          // Set the name of the TransportAction so it can fetch the appropriate material
          // properties.
          params.set<std::string>("transport_system") = name();
          // Group index and the number of groups are required to fetch the
          // scattering cross-section moments.
          params.set<unsigned int>("group_index") = g;
          params.set<unsigned int>("num_groups") = _num_groups;
          // Ordinate index is required to fetch the particle direction.
          params.set<unsigned int>("ordinate_index") = n;
          // Maximum scattering anisotropy.
          params.set<unsigned int>("max_anisotropy") = _max_eval_anisotropy;

          // Apply the parameters for the quadrature rule.
          applyQuadratureParameters(params);

          // Copy all of the group flux moment names into the variable
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
            params.set<std::vector<SubdomainName>>("block") =
                getParam<std::vector<SubdomainName>>("block");
          }

          _problem->addKernel("SAAFMomentScattering", "SAAFMomentScattering_" + var_name, params);
          debugOutput("      - Adding kernel SAAFMomentScattering for the variable " + var_name +
                      ".");
        } // SAAFMomentScattering
      }
    }
  }
}

void
TransportAction::addSAAFDiracKernels(const std::string & var_name, unsigned int g, unsigned int n)
{
  // Add SAAFPointSource.
  if (!_is_eigen && !_using_uncollided)
  {
    for (unsigned int i = 0u; i < _point_source_moments.size(); ++i)
    {
      auto params = _factory.getValidParams("SAAFPointSource");
      params.set<NonlinearVariableName>("variable") = var_name;
      // Set the name of the TransportAction so it can fetch the appropriate material properties.
      params.set<std::string>("transport_system") = name();
      params.set<MooseEnum>("point_not_found_behavior") =
          MooseEnum("ERROR WARNING IGNORE", "WARNING");
      // Group index is required to fetch the group particle removal cross-section
      // for stabilization.
      params.set<unsigned int>("group_index") = g;
      // Number of groups
      params.set<unsigned int>("num_groups") = _num_groups;
      // Ordinate index is required to fetch the particle direction.
      params.set<unsigned int>("ordinate_index") = n;

      params.set<Point>("point") = _point_source_locations[i];
      params.set<std::vector<Real>>("group_source") = _point_source_moments[i];
      params.set<unsigned int>("source_anisotropy") = _point_source_anisotropy[i];

      // Apply the parameters for the quadrature rule.
      applyQuadratureParameters(params);

      if (isParamValid("block"))
      {
        params.set<std::vector<SubdomainName>>("block") =
            getParam<std::vector<SubdomainName>>("block");
      }

      _problem->addDiracKernel("SAAFPointSource", "SAAFPointSource_" + var_name, params);
      debugOutput("      - Adding Dirac kernel SAAFPointSource for the "
                  "variable " +
                  var_name + ".");
    }
  } // SAAFPointSource
}
//------------------------------------------------------------------------------

//------------------------------------------------------------------------------
// Functions to add MOOSE objects for the diffusion approximation scheme.
//------------------------------------------------------------------------------
void
TransportAction::addDiffusionBCs(const std::string & var_name)
{
  // Add DiffusionRobinBC for vacuum boundary conditions.
  if (_vacuum_side_sets.size() > 0u)
  {
    auto params = _factory.getValidParams("DiffusionRobinBC");
    params.set<NonlinearVariableName>("variable") = var_name;
    params.set<Real>("incoming_partial_current") = 0.0;
    params.set<Real>("boundary_transport_correction") = 1.0;

    params.set<std::vector<BoundaryName>>("boundary") = _vacuum_side_sets;

    _problem->addBoundaryCondition("DiffusionRobinBC", "DiffusionRobinBC_" + var_name, params);
    debugOutput("      - Adding BC DiffusionRobinBC for the variable " + var_name + ".");
  } // DiffusionRobinBC

  // Add DiffusionRobinBC for surface source boundary conditions.
  if (_source_side_sets.size() > 0u)
  {
    mooseError("The diffusion approximation scheme does not support surface source boundary "
               "conditions.");

    auto params = _factory.getValidParams("DiffusionRobinBC");
    params.set<NonlinearVariableName>("variable") = var_name;
    // params.set<Real>("incoming_partial_current") = 0.0; // TODO: Implement me.
    params.set<Real>("boundary_transport_correction") = 1.0;

    params.set<std::vector<BoundaryName>>("boundary") = _source_side_sets;

    _problem->addBoundaryCondition("DiffusionRobinBC", "DiffusionRobinBC_" + var_name, params);
    debugOutput("      - Adding BC DiffusionRobinBC for the variable " + var_name + ".");
  } // DiffusionRobinBC

  // Add DiffusionNeumannBC for reflective boundaries.
  if (_reflective_side_sets.size() > 0u)
  {
    auto params = _factory.getValidParams("DiffusionNeumannBC");
    params.set<NonlinearVariableName>("variable") = var_name;
    params.set<Real>("net_partial_current") = 0.0;

    params.set<std::vector<BoundaryName>>("boundary") = _reflective_side_sets;

    _problem->addBoundaryCondition("DiffusionNeumannBC", "DiffusionNeumannBC_" + var_name, params);
    debugOutput("Adding BC DiffusionNeumannBC for the variable " + var_name + ".");
  } // DiffusionNeumannBC
}
void
TransportAction::addDiffusionICs(const std::string & var_name, unsigned int g)
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

        const auto & const_ic = getParam<std::vector<Real>>("constant_ic");
        if (const_ic.size() == 1)
          params.set<Real>("value") = const_ic[0];
        else if (const_ic.size() == _num_groups)
          params.set<Real>("value") = const_ic[g];
        else
        {
          mooseError("Size of 'constant_ic' does not match the declared number "
                     "of particle groups.");
        }

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
TransportAction::addDiffusionKernels(const std::string & var_name, unsigned int g)
{
  // Add ADParticleTimeDerivative.
  if (_exec_type == ExecutionType::Transient)
  {
    auto params = _factory.getValidParams("ADParticleTimeDerivative");
    params.set<NonlinearVariableName>("variable") = var_name;
    // Set the name of the TransportAction so it can fetch the appropriate material properties.
    params.set<std::string>("transport_system") = name();
    // Group index is required to fetch the group particle velocity.
    params.set<unsigned int>("group_index") = g;

    if (isParamValid("block"))
    {
      params.set<std::vector<SubdomainName>>("block") =
          getParam<std::vector<SubdomainName>>("block");
    }

    _problem->addKernel("ADParticleTimeDerivative", "ADParticleTimeDerivative_" + var_name, params);
    debugOutput("      - Adding kernel ADParticleTimeDerivative for the variable " + var_name +
                ".");
  } // ADParticleTimeDerivative

  // Add DiffusionApprox.
  {
    auto params = _factory.getValidParams("DiffusionApprox");
    params.set<NonlinearVariableName>("variable") = var_name;
    // Set the name of the TransportAction so it can fetch the appropriate material properties.
    params.set<std::string>("transport_system") = name();
    // Group index is required to fetch the group particle diffusion coefficient.
    params.set<unsigned int>("group_index") = g;

    if (isParamValid("block"))
    {
      params.set<std::vector<SubdomainName>>("block") =
          getParam<std::vector<SubdomainName>>("block");
    }

    _problem->addKernel("DiffusionApprox", "DiffusionApprox_" + var_name, params);
    debugOutput("      - Adding kernel DiffusionApprox for the variable " + var_name + ".");
  } // DiffusionApprox

  // Add DiffusionVolumeSource
  if (_volumetric_source_blocks.size() > 0u && !_is_eigen)
  {
    for (unsigned int i = 0u; i < _volumetric_source_blocks.size(); ++i)
    {
      auto params = _factory.getValidParams("DiffusionVolumeSource");
      params.set<NonlinearVariableName>("variable") = var_name;
      // Group index is required to fetch the group particle removal cross-section
      // for stabilization.
      params.set<unsigned int>("group_index") = g;
      // Ordinate index is required to fetch the particle direction.
      params.set<unsigned int>("num_groups") = _num_groups;
      params.set<std::vector<Real>>("group_source") = _volumetric_source_moments[i];

      params.set<std::vector<SubdomainName>>("block").emplace_back(_volumetric_source_blocks[i]);

      _problem->addKernel("DiffusionVolumeSource",
                          "DiffusionVolumeSource_" + var_name + "_" + _volumetric_source_blocks[i],
                          params);
      debugOutput("      - Adding kernel DiffusionVolumeSource for the variable " + var_name + ".");
    }
  } // DiffusionVolumeSource

  // Add DiffusionMaterialSource.
  if (!_is_eigen)
  {
    auto params = _factory.getValidParams("DiffusionMaterialSource");
    params.set<NonlinearVariableName>("variable") = var_name;
    // Set the name of the TransportAction so it can fetch the appropriate material properties.
    params.set<std::string>("transport_system") = name();
    // Group index and the number of groups are required to fetch the
    // source moments for the proper spectral energy group.
    params.set<unsigned int>("group_index") = g;
    params.set<unsigned int>("num_groups") = _num_groups;

    if (isParamValid("block"))
    {
      params.set<std::vector<SubdomainName>>("block") =
          getParam<std::vector<SubdomainName>>("block");
    }

    _problem->addKernel("DiffusionMaterialSource", "DiffusionMaterialSource_" + var_name, params);
    debugOutput("      - Adding kernel DiffusionMaterialSource for the variable " + var_name + ".");
  } // DiffusionMaterialSource

  // Only add fission kernels if debug doesn't disable them AND this transport system represents a
  // neutron field.
  if (!getParam<bool>("debug_disable_fission") && _particle == Particletype::Neutron)
  {
    // Add DiffusionFission.
    {
      auto params = _factory.getValidParams("DiffusionFission");
      params.set<NonlinearVariableName>("variable") = var_name;
      // Set the name of the TransportAction so it can fetch the appropriate material
      // properties.
      params.set<std::string>("transport_system") = name();
      // Group index and the number of groups are required to fetch the
      // scattering cross-section moments.
      params.set<unsigned int>("group_index") = g;
      params.set<unsigned int>("num_groups") = _num_groups;

      // Copy all of the scalar flux names into the variable parameter.
      auto & scalar_flux_names = params.set<std::vector<VariableName>>("group_scalar_fluxes");
      for (unsigned int g_prime = 0; g_prime < _num_groups; ++g_prime)
        scalar_flux_names.emplace_back(_group_flux_moments[g_prime][0u]);

      if (isParamValid("block"))
      {
        params.set<std::vector<SubdomainName>>("block") =
            getParam<std::vector<SubdomainName>>("block");
      }

      // For eigenvalues.
      if (_is_eigen)
        params.set<std::vector<TagName>>("extra_vector_tags").emplace_back("eigen");

      _problem->addKernel("DiffusionFission", "DiffusionFission_" + var_name, params);
      debugOutput("      - Adding kernel DiffusionFission for the variable " + var_name + ".");
    } // DiffusionFission
  }

  // Add DiffusionRemoval.
  {
    auto params = _factory.getValidParams("DiffusionRemoval");
    params.set<NonlinearVariableName>("variable") = var_name;
    // Set the name of the TransportAction so it can fetch the appropriate material properties.
    params.set<std::string>("transport_system") = name();
    // Group index is required to fetch the group particle diffusion coefficient.
    params.set<unsigned int>("group_index") = g;

    if (isParamValid("block"))
    {
      params.set<std::vector<SubdomainName>>("block") =
          getParam<std::vector<SubdomainName>>("block");
    }

    _problem->addKernel("DiffusionRemoval", "DiffusionRemoval_" + var_name, params);
    debugOutput("      - Adding kernel DiffusionRemoval for the variable " + var_name + ".");
  } // DiffusionRemoval

  // Add DiffusionScattering.
  if (_num_groups > 1u)
  {
    auto params = _factory.getValidParams("DiffusionScattering");
    params.set<NonlinearVariableName>("variable") = var_name;
    // Set the name of the TransportAction so it can fetch the appropriate material properties.
    params.set<std::string>("transport_system") = name();
    // Group index and the number of groups are required to fetch the
    // scattering cross-section moments.
    params.set<unsigned int>("group_index") = g;
    params.set<unsigned int>("num_groups") = _num_groups;

    // Copy all of the group flux ordinate names into the variable
    // parameter.
    auto & scalar_flux_names = params.set<std::vector<VariableName>>("group_scalar_fluxes");
    for (unsigned int g_prime = 0; g_prime < _num_groups; ++g_prime)
      scalar_flux_names.emplace_back(_group_flux_moments[g_prime][0u]);

    if (isParamValid("block"))
    {
      params.set<std::vector<SubdomainName>>("block") =
          getParam<std::vector<SubdomainName>>("block");
    }

    _problem->addKernel("DiffusionScattering", "DiffusionScattering_" + var_name, params);
    debugOutput("      - Adding kernel DiffusionScattering for the variable " + var_name + ".");
  } // DiffusionScattering
}
void
TransportAction::addDiffusionDiracKernels(const std::string & var_name, unsigned int g)
{
  // Add DiffusionIsoPointSource.
  if (!_is_eigen)
  {
    for (unsigned int i = 0u; i < _point_source_moments.size(); ++i)
    {
      auto params = _factory.getValidParams("DiffusionIsoPointSource");
      params.set<NonlinearVariableName>("variable") = var_name;
      params.set<MooseEnum>("point_not_found_behavior") =
          MooseEnum("ERROR WARNING IGNORE", "WARNING");
      // Group index is required to fetch the group particle removal cross-section
      // for stabilization.
      params.set<unsigned int>("group_index") = g;
      // Number of groups
      params.set<unsigned int>("num_groups") = _num_groups;

      params.set<Point>("point") = _point_source_locations[i];
      params.set<std::vector<Real>>("group_source") = _point_source_moments[i];
      params.set<unsigned int>("source_anisotropy") = _point_source_anisotropy[i];

      if (isParamValid("block"))
      {
        params.set<std::vector<SubdomainName>>("block") =
            getParam<std::vector<SubdomainName>>("block");
      }

      _problem->addDiracKernel(
          "DiffusionIsoPointSource", "DiffusionIsoPointSource_" + var_name, params);
      debugOutput("      - Adding Dirac kernel DiffusionIsoPointSource for the "
                  "variable " +
                  var_name + ".");
    }
  } // DiffusionIsoPointSource
}
//------------------------------------------------------------------------------
