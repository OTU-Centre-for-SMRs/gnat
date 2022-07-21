#include "NeutronTransportAction.h"

#include "Factory.h"
#include "Parser.h"
#include "NonlinearSystemBase.h"
#include "FEProblemBase.h"
#include "DGKernelBase.h"

#include "Conversion.h"
#include "MooseTypes.h"
#include "FEProblem.h"

#include "AddVariableAction.h"

#include "libmesh/string_to_enum.h"
#include "libmesh/fe_type.h"

registerMooseAction("GnatApp", NeutronTransportAction, "add_variable"); //
registerMooseAction("GnatApp", NeutronTransportAction, "add_kernel"); //
registerMooseAction("GnatApp", NeutronTransportAction, "add_dg_kernel"); //
registerMooseAction("GnatApp", NeutronTransportAction, "add_dirac_kernel"); //
registerMooseAction("GnatApp", NeutronTransportAction, "add_bc"); //
registerMooseAction("GnatApp", NeutronTransportAction, "add_ic"); //
registerMooseAction("GnatApp", NeutronTransportAction, "add_aux_variable"); //
registerMooseAction("GnatApp", NeutronTransportAction, "add_aux_kernel"); //

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
  params.addParam<std::vector<SubdomainName>>("block",
                                              "The list of blocks (ids or "
                                              "names) that this variable will "
                                              "be applied.");

  //----------------------------------------------------------------------------
  // Basic parameters for the neutron transport simulation.
  params.addRequiredRangeCheckedParam<unsigned int>("num_groups",
                                                    "num_groups > 0",
                                                    "The number of spectral "
                                                    "energy groups in the "
                                                    "problem.");
  params.addRequiredParam<MooseEnum>("execution_type", MooseEnum("steady transient"),
                                     "The method of execution for the problem. "
                                     "Options are steady-state source driven "
                                     "problems and transient source problems.");
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
  params.addParam<unsigned int>("max_anisotropy", 0,
                                "The maximum degree of anisotropy to evaluate. "
                                "Defaults to 0 for isotropic scattering.");
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

  //----------------------------------------------------------------------------
  // Parameters for debugging.
  params.addParam<bool>("debug_disable_scattering", false,
                        "Debug option to disable scattering evaluation.");
  params.addParam<bool>("debug_disable_source_iteration", true,
                        "Debug option to disable source iteration.");
  params.addParam<Real>("debug_steady_state_ic", 1.0,
                        "Debug initial guess for Newton's method.");

  return params;
}

NeutronTransportAction::NeutronTransportAction(const InputParameters & params)
  : Action(params)
  , _num_groups(getParam<unsigned int>("num_groups"))
  , _exec_type(getParam<MooseEnum>("execution_type").getEnum<ExecutionType>())
  , _n_l(getParam<unsigned int>("n_polar"))
  , _n_c(getParam<unsigned int>("n_azimuthal"))
  , _max_eval_anisotropy(getParam<unsigned int>("max_anisotropy"))
  , _angular_flux_name(getParam<std::string>("angular_flux_names"))
  , _flux_moment_name(getParam<std::string>("flux_moment_names"))
  , _vacuum_side_sets(getParam<std::vector<BoundaryName>>("vacuum_boundaries"))
  , _source_side_sets(getParam<std::vector<BoundaryName>>("source_boundaries"))
  , _reflective_side_sets(getParam<std::vector<BoundaryName>>("reflective_boundaries"))
  , _var_init(false)
{ }

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

void
NeutronTransportAction::addVariable(const std::string & var_name)
{
  auto fe_type = AddVariableAction::feType(_pars);
  auto type = AddVariableAction::variableType(fe_type, false, false);
  auto var_params = _factory.getValidParams(type);

  var_params.applySpecificParameters(_pars, {"family", "order"});
  var_params.set<std::vector<Real>>("scaling") = std::vector<Real>(getParam<Real>("scaling"));

  if (isParamValid("block"))
  {
    var_params.set<std::vector<SubdomainName>>("block")
      = getParam<std::vector<SubdomainName>>("block");
  }

  _problem->addVariable(type, var_name, var_params);
}

void
NeutronTransportAction::addKernels(const std::string & var_name, unsigned int g,
                                   unsigned int n)
{
  // Add ADNeutronTimeDerivative.
  if (_exec_type == ExecutionType::Transient)
  {
    auto params = _factory.getValidParams("ADNeutronTimeDerivative");
    params.set<NonlinearVariableName>("variable") = var_name;
    // Group index is required to fetch the group neutron velocity.
    params.set<unsigned int>("group_index") = g;

    if (isParamValid("block"))
    {
      params.set<std::vector<SubdomainName>>("block")
        = getParam<std::vector<SubdomainName>>("block");
    }

    _problem->addKernel("NtTimeDerivative",
                        "ADNeutronTimeDerivative_" + var_name,
                        params);
  } // ADNeutronTimeDerivative

  // Add ADNeutronStreaming.
  {
    auto params = _factory.getValidParams("ADNeutronStreaming");
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

    _problem->addKernel("ADNeutronStreaming",
                        "ADNeutronStreaming_" + var_name,
                        params);
  } // ADNeutronStreaming

  // Add ADNeutronRemoval.
  {
    auto params = _factory.getValidParams("ADNeutronRemoval");
    params.set<NonlinearVariableName>("variable") = var_name;
    // Group index is required to fetch the group neutron removal
    // cross-section.
    params.set<unsigned int>("group_index") = g;

    if (isParamValid("block"))
    {
      params.set<std::vector<SubdomainName>>("block")
        = getParam<std::vector<SubdomainName>>("block");
    }

    _problem->addKernel("ADNeutronRemoval",
                        "ADNeutronRemoval_" + var_name,
                        params);
  } // ADNeutronRemoval

  // Add ADNeutronMaterialSource.
  {
    auto params = _factory.getValidParams("ADNeutronMaterialSource");
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

    _problem->addKernel("ADNeutronMaterialSource",
                        "ADNeutronMaterialSource_" + var_name,
                        params);
  } // ADNeutronMaterialSource

  // Only add scattering kernels if debug doesn't disable them.
  if (!getParam<bool>("debug_disable_scattering"))
  {
    // Debug option to disable the source iteration solver.
    if (!getParam<bool>("debug_disable_source_iteration"))
    {
      // Compute the scattering evaluation with source iteration.
      mooseError("Scattering iteration currently is not supported and is a work "
                 "in progress.");

      // Add ADNeutronGToGScattering. This kernel is only enabled if there's
      // more than one group.
      if (_num_groups > 1u)
      {
        auto params = _factory.getValidParams("ADNeutronGToGScattering");
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

        _problem->addKernel("ADNeutronGToGScattering",
                            "ADNeutronGToGScattering_" + var_name,
                            params);
      } // ADNeutronGToGScattering

      // Add ADNeutronInGroupScattering.
      {
        auto params = _factory.getValidParams("ADNeutronInGroupScattering");
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

        _problem->addKernel("ADNeutronInGroupScattering",
                            "ADNeutronInGroupScattering_" + var_name,
                            params);
      } // ADNeutronInGroupScattering
    }
    else
    {
      // Computes the scattering evaluation without source iteration.
      // Add ADNeutronScattering.
      {
        auto params = _factory.getValidParams("ADNeutronScattering");
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

        _problem->addKernel("ADNeutronScattering",
                            "ADNeutronScattering_" + var_name,
                            params);
      } // ADNeutronScattering
    }
  }
}

void
NeutronTransportAction::addDGKernels(const std::string & var_name,
                                     unsigned int g, unsigned int n)
{
  // Add ADDGNeutronStreamingUpwind.
  {
    auto params = _factory.getValidParams("ADDGNeutronStreamingUpwind");
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

    _problem->addDGKernel("ADDGNeutronStreamingUpwind",
                          "ADDGNeutronStreamingUpwind_" + var_name,
                          params);
  } // ADDGNeutronStreamingUpwind
}

void
NeutronTransportAction::addBCs(const std::string & var_name, unsigned int g,
                               unsigned int n)
{
  // Add ADNeutronVacuumBC.
  if (_vacuum_side_sets.size() > 0u)
  {
    InputParameters params = _factory.getValidParams("ADNeutronVacuumBC");
    params.set<NonlinearVariableName>("variable") = var_name;
    // Ordinate index is required to fetch the neutron direction.
    params.set<unsigned int>("ordinate_index") = n;

    // Apply the parameters for the quadrature rule.
    applyQuadratureParameters(params);

    params.set<std::vector<BoundaryName>>("boundary") = _vacuum_side_sets;

    _problem->addBoundaryCondition("ADNeutronVacuumBC",
                                   "ADNeutronVacuumBC_" + var_name,
                                   params);
  } // ADNeutronVacuumBC

  // Add ADNeutronMatSourceBC.
  if (_source_side_sets.size() > 0u)
  {
    InputParameters params = _factory.getValidParams("ADNeutronMatSourceBC");
    params.set<NonlinearVariableName>("variable") = var_name;
    // Ordinate index is required to fetch the neutron direction.
    params.set<unsigned int>("ordinate_index") = n;

    // Apply the parameters for the quadrature rule.
    applyQuadratureParameters(params);

    params.set<std::vector<BoundaryName>>("boundary") = _source_side_sets;

    _problem->addBoundaryCondition("ADNeutronMatSourceBC",
                                   "ADNeutronMatSourceBC" + var_name,
                                   params);
  } // ADNeutronMatSourceBC

  // Add ADNeutronReflectiveBC.
  // TODO: Complete this implementation.
  if (_reflective_side_sets.size() > 0u)
  {
    mooseError("Reflective boundary conditions are currently not supported.");
  } // ADNeutronReflectiveBC
}

// TODO: Initial conditions.
void
NeutronTransportAction::addICs(const std::string & var_name, unsigned int g,
                               unsigned int n)
{
  // Set a constant initial condition for now.
  auto params = _factory.getValidParams("ConstantIC");
  params.set<VariableName>("variable") = var_name;

  if (isParamValid("block"))
  {
    params.set<std::vector<SubdomainName>>("block")
      = getParam<std::vector<SubdomainName>>("block");
  }

  params.set<Real>("value") = getParam<Real>("debug_steady_state_ic");
  _problem->addInitialCondition("ConstantIC", "ADNeutronMatSourceBC_" + var_name, params);
}

void
NeutronTransportAction::addAuxVariables(const std::string & var_name)
{
  auto fe_type = AddVariableAction::feType(_pars);
  auto type_name = AddVariableAction::variableType(fe_type, false, false);
  auto params = _factory.getValidParams(type_name);
  params.set<MooseEnum>("order") = fe_type.order.get_order();
  params.set<MooseEnum>("family") = Moose::stringify(fe_type.family);

  if (isParamValid("block"))
  {
    params.set<std::vector<SubdomainName>>("block")
      = getParam<std::vector<SubdomainName>>("block");
  }

  _problem->addAuxVariable(type_name, var_name, params);
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
  } // NeutronFluxMoment
}

void
NeutronTransportAction::addDiracKernels()
{
  // Add IsotropicNeutronPointSource.
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

        auto params = _factory.getValidParams("IsotropicNeutronPointSource");
        params.set<NonlinearVariableName>("variable") = var_name;
        params.set<MooseEnum>("point_not_found_behavior")
          = MooseEnum("ERROR WARNING IGNORE", "WARNING");

        // Apply the parameters for the quadrature rule.
        applyQuadratureParameters(params);

        // Set the intensities and points for the group g.
        auto & points = params.set<std::vector<Point>>("points");
        auto & intensities = params.set<std::vector<Real>>("intensities");
        points.reserve(mapping.size());
        intensities.reserve(mapping.size());
        for (const auto & index : mapping)
        {
          points.emplace_back(source_locations[index]);
          intensities.emplace_back(source_intensities[index]);
        }

        if (isParamValid("block"))
        {
          params.set<std::vector<SubdomainName>>("block")
            = getParam<std::vector<SubdomainName>>("block");
        }

        _problem->addDiracKernel("IsotropicNeutronPointSource",
                                 "IsotropicNeutronPointSource_" + var_name,
                                 params);
      }
    }
  } // IsotropicNeutronPointSource
}

void
NeutronTransportAction::act()
{
  if (!_var_init)
  {
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
  }

  // Loop over all groups.
  for (unsigned int g = 0; g < _num_groups; ++g)
  {
    // Loop over all the flux ordinates.
    for (unsigned int n = 0; n < _num_flux_ordinates; ++n)
    {
      const auto & var_name = _group_angular_fluxes[g][n];

      // Add a non-linear variable.
      if (_current_task == "add_variable")
        addVariable(var_name);

      // Add kernels.
      if (_current_task == "add_kernel")
        addKernels(var_name, g, n);

      // Add DG kernels.
      if (_current_task == "add_dg_kernel")
        addDGKernels(var_name, g, n);

      // Add boundary conditions.
      if (_current_task == "add_bc")
        addBCs(var_name, g, n);

      // Add initial conditions.
      if (_current_task == "add_ic")
        addICs(var_name, g, n);
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
          addAuxVariables(var_name);

        // Add auxkernels.
        if (_current_task == "add_aux_kernel")
          addAuxKernels(var_name, g, l, m);

        moment_index++;
      }
    }
  }

  // Add Dirac kernels (isotropic point sources). Anisotropic point sources are
  // handled separately.
  if (_current_task == "add_dirac_kernel")
    addDiracKernels();
}

void
NeutronTransportAction::addRelationshipManagers(Moose::RelationshipManagerType input_rm_type)
{
  auto dg_kernel_params = DGKernelBase::validParams();
  addRelationshipManagers(input_rm_type, dg_kernel_params);
}
