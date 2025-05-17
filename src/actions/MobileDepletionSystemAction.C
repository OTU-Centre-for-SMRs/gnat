#include "MobileDepletionSystemAction.h"

#include <filesystem>
#include <algorithm>

#include "AddVariableAction.h"
#include "AddOutputAction.h"
#include "TransportAction.h"
#include "InputParameterWarehouse.h"

#include "WCNSFVFlowPhysics.h"

#include "DepletionLibraryAction.h"

registerMooseAction("GnatApp", MobileDepletionSystemAction, "add_variable");
registerMooseAction("GnatApp", MobileDepletionSystemAction, "add_ic");
registerMooseAction("GnatApp", MobileDepletionSystemAction, "add_material");

// Add auxvariables and kernels for nuclide mass/number densities since we solve for the mass
// fraction as the primal unknown.
registerMooseAction("GnatApp", MobileDepletionSystemAction, "add_aux_variable");
registerMooseAction("GnatApp", MobileDepletionSystemAction, "add_aux_kernel");

// For the SUPG finite element scheme.
registerMooseAction("GnatApp", MobileDepletionSystemAction, "add_kernel");
registerMooseAction("GnatApp", MobileDepletionSystemAction, "add_bc");

// For the finite volume scheme.
registerMooseAction("GnatApp", MobileDepletionSystemAction, "add_fv_kernel");
registerMooseAction("GnatApp", MobileDepletionSystemAction, "add_fv_bc");

InputParameters
MobileDepletionSystemAction::validParams()
{
  auto params = GnatBaseAction::validParams();
  params.addClassDescription(
      "This action adds variables which represents the density of all mobile nuclides "
      "contained in the provided cross-section files, alongside the kernels required to describe "
      "the time evolution of those nuclides.");

  //----------------------------------------------------------------------------
  // Parameters required for mass transport.
  params.addRequiredParam<MooseEnum>(
      "scheme",
      MooseEnum("supg_fe fv"),
      "The discretization and stabilization scheme that the nuclide system should use.");

  // The coupled TransportSystem.
  params.addParam<std::string>(
      "transport_system", "", "Name of the transport system which will provide scalar fluxes.");

  //----------------------------------------------------------------------------
  // Parameters for finite element simulations.
  params.addParam<MooseEnum>("family",
                             AddVariableAction::getNonlinearVariableFamilies(),
                             "Specifies the family of FE shape functions to "
                             "use for this variable.");
  params.addParam<MooseEnum>("order",
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
  // Parameters for finite volume simulations.
  params.addParam<bool>("using_moose_ns_fv",
                        false,
                        "Whether the simulation should assume that the MOOSE Navier-Stokes module "
                        "is being used for flow capabilities. This allows us to use more "
                        "sophisticated stabilization schemes.");
  params.addParam<MooseEnum>(
      "fv_face_interpolation",
      MooseEnum("average skewness-corrected", "average"),
      "The numerical scheme to interpolate the nuclide scalar field variables to the "
      "face (separate from the advected quantity interpolation).");
  params.addParam<bool>(
      "fv_two_term_boundary_expansion",
      false,
      "Whether the simulation should use a two-term Taylor series expansion on the boundaries.");
  params.addParam<MooseEnum>(
      "fv_adv_interpolation",
      MooseEnum("average upwind skewness-corrected min_mod vanLeer", "average"),
      "The numerical scheme to use for interpolating nuclide scalar fields, as an advected "
      "quantity, to the face.");

  //----------------------------------------------------------------------------
  // Parameters required for advection.
  params.addParam<MooseFunctorName>("u", "The functor name for the x-component of the velocity.");
  params.addParam<MooseFunctorName>("v", "The functor name for the y-component of the velocity.");
  params.addParam<MooseFunctorName>("w", "The functor name for the z-component of the velocity.");

  //----------------------------------------------------------------------------
  // Diffusion coefficient parameters.
  params.addParam<MooseEnum>("turbulence_handling",
                             MooseEnum("none mixing-length"),
                             "What type of turbulent diffusion to use.");
  params.addParam<MooseFunctorName>("density", "The density of the bulk fluid ($g/cm^{3}$");
  params.addParam<MooseFunctorName>("temperature", "The temperature of the bulk fluid ($K$).");
  params.addParam<MooseFunctorName>("dynamic_viscosity",
                                    "The dynamic viscosity of the bulk fluid ($g/(cm s)$).");

  //----------------------------------------------------------------------------
  // Turbulence modelling parameters.
  params.addParam<Real>("schmidt_number",
                        0.7,
                        "The turbulent Schmidt number that relates the turbulent scalar "
                        "diffusivity to the turbulent momentum diffusivity.");
  params.addParam<MooseFunctorName>("mixing_length", "The name of the mixing length functor.");

  //----------------------------------------------------------------------------
  // Boundary condition parameters.
  params.addParam<std::vector<BoundaryName>>(
      "inlet_boundaries",
      std::vector<BoundaryName>(),
      "The names of boundaries which act as inflow boundaries. GNAT applies a fixed nuclide mass "
      "fraction at each of these boundaries equal to the initial mass fractions specified in "
      "'element_atom_fractions' and 'extra_nuclide_atom_fractions'. This functionality can be "
      "overriden by specifying inlet mass fractions in 'inlet_mass_fractions'.");
  params.addParam<std::vector<std::vector<Real>>>(
      "inlet_atom_fractions",
      std::vector<std::vector<Real>>(),
      "A double vector of inlet atom fractions. The outer vector must match the ordering of "
      "'inlet_boundaries'. The inner vector must must be arranged so elements go first; extra "
      "nuclides go second. The ordering of elements and nuclides must match the ordering of "
      "'elements' and 'extra_nuclides'.");
  params.addParam<std::vector<BoundaryName>>(
      "outlet_boundaries",
      std::vector<BoundaryName>(),
      "The names of boundaries which act as outflow boundaries. If 'using_moose_ns_fv' is set to "
      "'true' this parameter will be ignored.");

  //----------------------------------------------------------------------------
  // Initial conditions and elemental composition parameters.
  params.addParam<std::vector<std::string>>(
      "elements",
      std::vector<std::string>(),
      "The elemental composition of the fluid mixture. It is assumed that all elements added with "
      "this syntax are composed of a mixture of nuclides at their natural abundances.");
  params.addParam<std::vector<Real>>(
      "element_atom_fractions",
      std::vector<Real>(),
      "The atom fractions of each element in 'elements'. Must "
      "follow the same ordering as 'elements'. The sum of all atom fractions (multiplied by the "
      "natural) provided must equal one.");
  params.addParam<std::vector<std::string>>(
      "extra_nuclides",
      std::vector<std::string>(),
      "Extra nuclides to add to the composition of the fluid mixture.");
  params.addParam<std::vector<Real>>("extra_nuclide_atom_fractions",
                                     std::vector<Real>(),
                                     "The atom fractions of each nuclide in 'extra_nuclides'. Must "
                                     "follow the same ordering as 'extra_nuclides'. The sum of all "
                                     "atom fractions provided must equal one.");
  params.addParam<bool>("compute_number_density",
                        true,
                        "Whether a number density should be computed or not. Setting this flag to "
                        "false results in a mass density instead.");
  params.addParam<bool>("output_nuclide_mass_fractions",
                        false,
                        "Whether the mass fractions of nuclides should be output alongside their "
                        "respective densities.");

  //----------------------------------------------------------------------------
  // Parameters for radiological consequence assessment simulations; coupled radiation transport
  // with mobile depletion where the depleted radionuclides act as photon/neutron decay sources.
  params.addParam<bool>("add_photon_sources",
                        false,
                        "If group-wise radionuclide photon source terms should be added.");
  params.addParam<std::string>(
      "photon_source_prefix",
      "photon_source",
      "A prefix for the auxvariables which stores the group-wise photon particle source.");
  params.addParam<bool>("add_neutron_sources",
                        false,
                        "If group-wise radionuclide neutron source terms should be added.");
  params.addParam<std::string>(
      "neutron_source_prefix",
      "neutron_source",
      "A prefix for the auxvariables which stores the group-wise neutron particle source.");

  params.addParam<std::vector<Real>>(
      "photon_group_boundaries",
      std::vector<Real>{2.0e7, 0.0},
      "The photon group structure in descending order of energy (including 0.0 eV)");
  params.addParam<std::vector<Real>>(
      "neutron_group_boundaries",
      std::vector<Real>{2.0e7, 0.0},
      "The neutron group structure in descending order of energy (including 0.0 eV)");

  //----------------------------------------------------------------------------
  // Debug parameters.
  params.addParam<std::vector<std::string>>(
      "debug_filter_nuclides",
      std::vector<std::string>(),
      "A list of nuclides that should be included in the mobile depletion analysis. ");

  return params;
}

MobileDepletionSystemAction::MobileDepletionSystemAction(const InputParameters & parameters)
  : GnatBaseAction(parameters),
    _mesh_dims(0u),
    _using_moose_ns(getParam<bool>("using_moose_ns_fv")),
    _scheme(getParam<MooseEnum>("scheme").getEnum<NuclideScheme>()),
    _coupled_ns_fv(nullptr),
    _coupled_depletion_lib(nullptr),
    _transport_system(getParam<std::string>("transport_system")),
    _has_transport_system(true),
    _num_groups(0u),
    _first_action(true),
    _elements(getParam<std::vector<std::string>>("elements")),
    _element_atom_fractions(getParam<std::vector<Real>>("element_atom_fractions")),
    _extra_nuclides(getParam<std::vector<std::string>>("extra_nuclides")),
    _extra_nuclide_atom_fractions(getParam<std::vector<Real>>("extra_nuclide_atom_fractions")),
    _inlet_boundaries(getParam<std::vector<BoundaryName>>("inlet_boundaries")),
    _inlet_atom_fractions(getParam<std::vector<std::vector<Real>>>("inlet_atom_fractions")),
    _outlet_boundaries(getParam<std::vector<BoundaryName>>("outlet_boundaries")),
    _photon_source_prefix(getParam<std::string>("photon_source_prefix")),
    _neutron_source_prefix(getParam<std::string>("neutron_source_prefix")),
    _photon_group_boundaries(getParam<std::vector<Real>>("photon_group_boundaries")),
    _neutron_group_boundaries(getParam<std::vector<Real>>("neutron_group_boundaries")),
    _debug_filter_nuclides(getParam<std::vector<std::string>>("debug_filter_nuclides"))
{
  if (_scheme == NuclideScheme::SUPGFE && !_pars.isParamSetByUser("family"))
    paramError("family", "'order' must be set for SUPG finite element dispersion simulations.");
  if (_scheme == NuclideScheme::SUPGFE && !_pars.isParamSetByUser("order"))
    paramError("order", "'family' must be set for SUPG finite element dispersion simulations.");

  // Error handle the user not providing the appropriate material properties.
  if (!_using_moose_ns)
  {
    if (!isParamValid("u"))
      paramError("u", "The x component of the velocity must be supplied using the 'u' parameter.");

    if (!isParamValid("density"))
      paramError("density",
                 "The material density must be specified if you are not using the MOOSE "
                 "Navier-Stokes finite volume system.");

    if (!isParamValid("temperature"))
      paramError("temperature",
                 "The material temperature must be specified if you are not using the MOOSE "
                 "Navier-Stokes finite volume system.");

    if (!isParamValid("dynamic_viscosity"))
      paramError("dynamic_viscosity",
                 "The material dynamic viscosity must be specified if you are not using the MOOSE "
                 "Navier-Stokes finite volume system.");

    if (!isParamValid("mixing_length") && getParam<MooseEnum>("turbulence_handling") != "none")
      paramError("mixing_length",
                 "A mixing length functor must be specified if you are not using the MOOSE "
                 "Navier-Stokes finite volume system.");
  }

  // Handle no elements or nuclides being provided.
  if (_elements.size() + _extra_nuclides.size() == 0u)
    mooseError("No elements or nuclides have been specified!");

  // Handle issues where there are element/nuclide mismatches with atom fractions.
  if (_elements.size() != _element_atom_fractions.size())
    mooseError("There is a mismatch between the number of provided elements and the number of "
               "provided element atom fractions.");
  if (_extra_nuclides.size() != _extra_nuclide_atom_fractions.size())
    mooseError("There is a mismatch between the number of provided nuclides and the number of "
               "provided nuclide atom fractions.");

  // Sanity check to make sure all user provided elements and nuclides are actually elements.
  for (const auto & element : _elements)
    if (!NuclearData::Nuclide::isElement(element))
      paramError("elements", "The element " + element + " is not a real element.");

  for (const auto & nuclide : _extra_nuclides)
    if (!NuclearData::Nuclide::isElement(nuclide))
      paramError("extra_nuclides",
                 "The nuclide " + nuclide + " is not a nuclide of a real element.");

  // Error handle malformed boundary conditions.
  if (_inlet_atom_fractions.size() > 0u)
  {
    if (_inlet_boundaries.size() != _inlet_atom_fractions.size())
      mooseError("The number of provided inlet boundaries does not match the number of provided "
                 "inlet mass fractions.");

    for (unsigned int i = 0u; i < _inlet_atom_fractions.size(); ++i)
    {
      if (_inlet_atom_fractions[i].size() != _elements.size() + _extra_nuclides.size())
        paramError("inlet_mass_fractions",
                   "The boundary '" + _inlet_boundaries[i] +
                       "' is missing element/nuclide atom fractions. " +
                       std::to_string(_inlet_atom_fractions[i].size()) +
                       " atom fractions have been provided, " +
                       std::to_string(_elements.size() + _extra_nuclides.size()) +
                       " are required.");
    }
  }
}

// Add elements using their natural isotopic abundance.
void
MobileDepletionSystemAction::addElementsAndNuclides()
{
  // First store everything in atom fractions.
  for (unsigned int i = 0u; i < _elements.size(); ++i)
  {
    auto abundances = NuclearData::Nuclide::getAbundances(_elements[i]);
    for (auto & [nuclide, abundance] : abundances)
    {
      if (_total_nuclide_list.count(nuclide) == 0u)
        _total_nuclide_list.emplace(nuclide, abundance * _element_atom_fractions[i]);
      else
        _console << COLOR_YELLOW << "The nuclide " << nuclide
                 << " has been provided multiple times. This iteration (element atom fraction of "
                 << _element_atom_fractions[i]
                 << ") will be ignored when computing the composition of the mixture.\n"
                 << COLOR_DEFAULT;
    }
  }

  for (unsigned int i = 0u; i < _extra_nuclides.size(); ++i)
  {
    if (_total_nuclide_list.count(_extra_nuclides[i]) == 0u)
      _total_nuclide_list.emplace(_extra_nuclides[i], _extra_nuclide_atom_fractions[i]);
    else
      _console << COLOR_YELLOW << "The nuclide " << _extra_nuclides[i]
               << " has been provided multiple times. This iteration (nuclide atom fraction of "
               << _extra_nuclide_atom_fractions[i]
               << ") will be ignored when computing the composition of the mixture.\n"
               << COLOR_DEFAULT;
  }

  // Now perform a normalization to compute weight fractions.
  {
    Real total_weight = 0.0;
    for (auto & [nuclide, fraction] : _total_nuclide_list)
    {
      fraction *= NuclearData::Nuclide::getAtomicMass(nuclide);
      total_weight += fraction;
    }
    if (!MooseUtils::absoluteFuzzyEqual(total_weight, 0.0))
      for (auto & [nuclide, fraction] : _total_nuclide_list)
        fraction /= total_weight;
  }

  // Sanity check to make sure everything adds up to 1.0. If not, warn the user and scale each
  // weight fraction to ensure the results are conservative.
  {
    Real sum = 0.0;
    for (auto & [nuclide, fraction] : _total_nuclide_list)
      sum += fraction;

    if (MooseUtils::absoluteFuzzyGreaterThan(sum, 1.0))
    {
      _console << COLOR_YELLOW << "The sum of all computed weight fractions (" << sum
               << ") is greater than 1. Each weight fraction will be divided by " << sum
               << " to remain conservative.\n"
               << COLOR_DEFAULT;
      for (auto & [nuclide, fraction] : _total_nuclide_list)
        fraction /= sum;
    }
  }
  _console << std::flush;

  // Repeat the same process for the boundary nuclide atom fractions.
  if (_inlet_atom_fractions.size() > 0u)
  {
    _boundary_mass_fractions.resize(_inlet_boundaries.size());

    // Loop over boundaries.
    for (unsigned int i = 0u; i < _inlet_boundaries.size(); ++i)
    {
      const auto & atom_fractions = _inlet_atom_fractions[i];

      // First loop over the element part of the array.
      for (unsigned int j = 0u; j < _elements.size(); ++j)
      {
        auto abundances = NuclearData::Nuclide::getAbundances(_elements[j]);
        for (const auto & [nuclide, abundance] : abundances)
        {
          if (_boundary_mass_fractions[i].count(nuclide) == 0u)
            _boundary_mass_fractions[i].emplace(nuclide, abundance * atom_fractions[j]);
        }
      }

      // Next, loop over the extra nuclide part.
      for (unsigned int j = 0u; j < _extra_nuclides.size(); ++j)
      {
        if (_boundary_mass_fractions[i].count(_extra_nuclides[j]) == 0u)
          _boundary_mass_fractions[i].emplace(_extra_nuclides[j],
                                              atom_fractions[_elements.size() + j]);
      }

      // Perform a normalization to compute weight fractions.
      {
        Real total_weight = 0.0;
        for (auto & [nuclide, fraction] : _boundary_mass_fractions[i])
        {
          fraction *= NuclearData::Nuclide::getAtomicMass(nuclide);
          total_weight += fraction;
        }
        if (!MooseUtils::absoluteFuzzyEqual(total_weight, 0.0))
          for (auto & [nuclide, fraction] : _boundary_mass_fractions[i])
            fraction /= total_weight;
      }

      // Sanity check to make sure everything adds up to 1.0. If not, warn the user and scale each
      // weight fraction to ensure the results are conservative.
      {
        Real sum = 0.0;
        for (auto & [nuclide, fraction] : _boundary_mass_fractions[i])
          sum += fraction;

        if (MooseUtils::absoluteFuzzyGreaterThan(sum, 1.0))
        {
          _console << COLOR_YELLOW << "The sum of all computed boundary weight fractions (" << sum
                   << ") is greater than 1. Each weight fraction will be divided by " << sum
                   << " to remain conservative.\n"
                   << COLOR_DEFAULT;
          for (auto & [nuclide, fraction] : _boundary_mass_fractions[i])
            fraction /= sum;
        }
      }
    }
  }
}

// Build the depletion chain and fetch depletion properties.
void
MobileDepletionSystemAction::buildDepletionSystem()
{
  std::vector<std::string> nuclides;
  nuclides.reserve(_total_nuclide_list.size());
  for (const auto & [nuclide, density] : _total_nuclide_list)
    nuclides.emplace_back(nuclide);

  // Fetch the entire depletion chain.
  _coupled_depletion_lib->getNuclidesFromInitial(nuclides);
  for (const auto & nuclide : nuclides)
  {
    if (_total_nuclide_list.count(nuclide) == 0u)
      _total_nuclide_list.emplace(nuclide, 0.0);
  }

  // Remove nuclides that do not exist in the depletion system.
  {
    std::vector<std::string> remove;
    for (const auto & [nuclide, weight] : _total_nuclide_list)
    {
      if (!_coupled_depletion_lib->hasNuclide(nuclide))
        remove.emplace_back(nuclide);
    }

    for (const auto & nuclide : remove)
      _total_nuclide_list.erase(nuclide);
  }

  // Prune the nuclide list if the debug filter is enabled.
  if (_debug_filter_nuclides.size() > 0u)
  {
    std::vector<std::string> to_remove;
    for (const auto & [nuclide, fraction] : _total_nuclide_list)
    {
      if (std::find(_debug_filter_nuclides.begin(), _debug_filter_nuclides.end(), nuclide) ==
          _debug_filter_nuclides.end())
        to_remove.emplace_back(nuclide);
    }
    for (const auto & remove_nuclide : to_remove)
      _total_nuclide_list.erase(remove_nuclide);
  }

  // Add these to the boundary condition mass fraction list.
  for (unsigned int i = 0u; i < _boundary_mass_fractions.size(); ++i)
  {
    for (const auto & [nuclide, fraction] : _total_nuclide_list)
    {
      if (_boundary_mass_fractions[i].count(nuclide) == 0.0)
        _boundary_mass_fractions[i].emplace(nuclide, fraction);
    }
  }

  // Sanity check to make sure the number of cross-section energy groups are the same.
  if (_num_groups != _coupled_depletion_lib->numGroups() && _has_transport_system)
    mooseError("The number of microscopic cross-section energy groups does not match the number of "
               "neutron energy groups.");
}

void
MobileDepletionSystemAction::applyIsotopeParameters(InputParameters & params, bool apply_density)
{
  if (_coupled_ns_fv)
  {
    if (apply_density)
      params.set<MooseFunctorName>("density") = _coupled_ns_fv->densityName();

    const auto & vel_names = _coupled_ns_fv->getVelocityNames();
    params.set<MooseFunctorName>("u") = vel_names[0];
    if (_mesh_dims >= 2u)
      params.set<MooseFunctorName>("v") = vel_names[1];
    if (_mesh_dims >= 3u)
      params.set<MooseFunctorName>("w") = vel_names[2];
  }
  else
  {
    if (apply_density)
      params.set<MooseFunctorName>("density") = getParam<MooseFunctorName>("density");

    params.set<MooseFunctorName>("u") = getParam<MooseFunctorName>("u");
    if (_mesh_dims >= 2u && isParamValid("v"))
      params.set<MooseFunctorName>("v") = getParam<MooseFunctorName>("v");
    else
      paramError(
          "v",
          "In >= 2D the y component of the velocity must be supplied using the 'v' parameter.");
    if (_mesh_dims >= 3u && isParamValid("w"))
      params.set<MooseFunctorName>("w") = getParam<MooseFunctorName>("w");
    else
      paramError("w",
                 "In 3D the z component of the velocity must be supplied using the 'w' parameter.");
  }
}

void
MobileDepletionSystemAction::addICs(const std::string & nuclide_var_name)
{
  const std::string frac_var_name = nuclide_var_name + "_mass_fraction";

  // Add ConstantIC.
  {
    auto params = _factory.getValidParams("ConstantIC");
    params.set<VariableName>("variable") = frac_var_name;

    if (isParamValid("block"))
    {
      params.set<std::vector<SubdomainName>>("block") =
          getParam<std::vector<SubdomainName>>("block");
    }

    params.set<Real>("value") = _total_nuclide_list[nuclide_var_name];

    _problem->addInitialCondition("ConstantIC", "ConstantIC_" + frac_var_name, params);
    debugOutput("Adding IC ConstantIC for the variable " + frac_var_name + ".");
  } // ConstantIC
}

void
MobileDepletionSystemAction::addMaterials(const std::string & nuclide_var_name)
{
  const std::string frac_var_name = nuclide_var_name + "_mass_fraction";

  // Add FunctorAutoNuclideMaterial.
  {
    auto params = _factory.getValidParams("FunctorAutoNuclideMaterial");
    params.set<NonlinearVariableName>("isotope_name") = frac_var_name;

    // Set parameters for the Brownian diffusion coefficient.
    params.set<MooseEnum>("scheme") = getParam<MooseEnum>("scheme");
    params.set<MooseEnum>("turbulence_handling") = getParam<MooseEnum>("turbulence_handling");
    params.set<Real>("radii") = 5.3e-9; // TODO: Actual radii of the particles.

    if (_coupled_ns_fv)
    {
      params.set<MooseFunctorName>("temperature") = _coupled_ns_fv->getFluidTemperatureName();

      params.set<MooseFunctorName>("dynamic_viscosity") =
          _coupled_ns_fv->getParam<MooseFunctorName>("dynamic_viscosity");
    }
    else
    {
      params.set<MooseFunctorName>("temperature") = getParam<MooseFunctorName>("temperature");
      params.set<MooseFunctorName>("dynamic_viscosity") =
          getParam<MooseFunctorName>("dynamic_viscosity");
    }

    // Set parameters for the turbulence modelling scheme.
    params.set<Real>("schmidt_number") = getParam<Real>("schmidt_number");
    if (_coupled_ns_fv && getParam<MooseEnum>("turbulence_handling") == "mixing-length")
      params.set<MooseFunctorName>("mixing_length") = "mixing_length";
    else
      params.set<MooseFunctorName>("mixing_length") = getParam<MooseFunctorName>("mixing_length");

    applyIsotopeParameters(params, false);

    _problem->addMaterial(
        "FunctorAutoNuclideMaterial", "FunctorAutoNuclideMaterial" + frac_var_name, params);
    debugOutput("    - Adding material FunctorAutoNuclideMaterial for the variable " +
                frac_var_name + ".");
  } // FunctorAutoNuclideMaterial
}

void
MobileDepletionSystemAction::addKernels(const std::string & nuclide_var_name)
{
  const std::string frac_var_name = nuclide_var_name + "_mass_fraction";

  const auto & nuclide_data = _coupled_depletion_lib->getNuclide(nuclide_var_name);
  const auto & nuclide_rxns = nuclide_data.getReactions();
  const auto & activation_parents = _coupled_depletion_lib->getActivationParents(nuclide_var_name);
  const auto & decay_parents = _coupled_depletion_lib->getDecayParents(nuclide_var_name);

  // Add ADMassFractionNuclideActivation.
  bool should_add_act_src = false;
  for (const auto & act_nuclide : activation_parents)
  {
    if (_total_nuclide_list.count(act_nuclide) > 0u)
    {
      should_add_act_src = true;
      break;
    }
  }
  if (should_add_act_src && _has_transport_system)
  {
    auto params = _factory.getValidParams("ADMassFractionNuclideActivation");
    params.set<NonlinearVariableName>("variable") = frac_var_name;
    // Set the number of groups.
    params.set<unsigned int>("num_groups") = _num_groups;

    // Set the activation cross-sections and parent mass fractions.
    auto & all_xs = params.set<std::vector<Real>>("group_activation");
    auto & fractions = params.set<std::vector<VariableName>>("isotope_mass_fractions");

    unsigned int i = 0u;
    for (const auto & act_nuclide : activation_parents)
    {
      if (_total_nuclide_list.count(act_nuclide) == 0u)
        continue;

      fractions.emplace_back(act_nuclide + "_mass_fraction");
      for (unsigned int g = 0u; g < _num_groups; ++g)
        all_xs.emplace_back(0.0);

      for (const auto & rxn : _coupled_depletion_lib->getNuclide(act_nuclide).getReactions())
      {
        if (rxn._cross_sections.size() > 0u)
        {
          for (unsigned int g = 0u; g < _num_groups; ++g)
            all_xs[i * _num_groups + g] += rxn._branching_factor * rxn._cross_sections[g];
        }
      }
      i++;
    }

    // Apply the scalar fluxes.
    auto & scalar_flux_names = params.set<std::vector<VariableName>>("group_scalar_fluxes");
    scalar_flux_names.reserve(_num_groups);
    for (unsigned int g = 0u; g < _num_groups; ++g)
      scalar_flux_names.emplace_back(_group_flux_moments[g]);

    // Apply common SUPG parameters.
    applyIsotopeParameters(params);

    if (isParamValid("block"))
    {
      params.set<std::vector<SubdomainName>>("block") =
          getParam<std::vector<SubdomainName>>("block");
    }

    _problem->addKernel("ADMassFractionNuclideActivation",
                        "ADMassFractionNuclideActivation" + frac_var_name,
                        params);
    debugOutput("    - Adding kernel ADMassFractionNuclideActivation for the variable " +
                frac_var_name + ".");
  } // ADMassFractionNuclideActivation

  // Add ADMassFractionNuclideAdvection.
  {
    auto params = _factory.getValidParams("ADMassFractionNuclideAdvection");
    params.set<NonlinearVariableName>("variable") = frac_var_name;

    // Apply common SUPG parameters.
    applyIsotopeParameters(params);

    if (isParamValid("block"))
    {
      params.set<std::vector<SubdomainName>>("block") =
          getParam<std::vector<SubdomainName>>("block");
    }

    _problem->addKernel(
        "ADMassFractionNuclideAdvection", "ADMassFractionNuclideAdvection" + frac_var_name, params);
    debugOutput("    - Adding kernel ADMassFractionNuclideAdvection for the variable " +
                frac_var_name + ".");
  }

  // Add ADMassFractionNuclideDecaySink.
  if (nuclide_data.decayConst() > 0.0)
  {
    auto params = _factory.getValidParams("ADMassFractionNuclideDecaySink");
    params.set<NonlinearVariableName>("variable") = frac_var_name;
    // Decay constant.
    params.set<Real>("decay_const") = nuclide_data.decayConst();

    // Apply common SUPG parameters.
    applyIsotopeParameters(params);

    if (isParamValid("block"))
    {
      params.set<std::vector<SubdomainName>>("block") =
          getParam<std::vector<SubdomainName>>("block");
    }

    _problem->addKernel(
        "ADMassFractionNuclideDecaySink", "ADMassFractionNuclideDecaySink" + frac_var_name, params);
    debugOutput("    - Adding kernel ADMassFractionNuclideDecaySink for the variable " +
                frac_var_name + ".");
  } // ADMassFractionNuclideDecaySink

  // Add ADMassFractionNuclideDecaySource.
  bool should_add_decay_src = false;
  for (const auto & dec_nuclide : decay_parents)
  {
    if (_total_nuclide_list.count(dec_nuclide) > 0u)
    {
      should_add_decay_src = true;
      break;
    }
  }
  if (should_add_decay_src)
  {
    auto params = _factory.getValidParams("ADMassFractionNuclideDecaySource");
    params.set<NonlinearVariableName>("variable") = frac_var_name;

    // Apply decay constants and branching factors.
    auto & decay_data = params.set<std::vector<Real>>("decay_constants");

    auto & fractions = params.set<std::vector<VariableName>>("isotope_mass_fractions");

    unsigned int i = 0u;
    for (const auto & dec_nuclide : decay_parents)
    {
      if (_total_nuclide_list.count(dec_nuclide) == 0u)
        continue;

      fractions.emplace_back(dec_nuclide + "_mass_fraction");
      decay_data.emplace_back(0.0);

      for (const auto & decay : _coupled_depletion_lib->getNuclide(dec_nuclide).getDecays())
      {
        if (decay._target == nuclide_var_name)
          decay_data[i] += decay._branching_factor;
      }
      decay_data[i] *= _coupled_depletion_lib->getNuclide(dec_nuclide).decayConst();
      i++;
    }

    // Apply common SUPG parameters.
    applyIsotopeParameters(params);

    if (isParamValid("block"))
    {
      params.set<std::vector<SubdomainName>>("block") =
          getParam<std::vector<SubdomainName>>("block");
    }

    _problem->addKernel("ADMassFractionNuclideDecaySource",
                        "ADMassFractionNuclideDecaySource" + frac_var_name,
                        params);
    debugOutput("    - Adding kernel ADMassFractionNuclideDecaySource for the variable " +
                frac_var_name + ".");
  } // ADMassFractionNuclideDecaySource

  // Add ADMassFractionNuclideDepletion.
  if (nuclide_rxns.size() > 0u && _has_transport_system)
  {
    auto params = _factory.getValidParams("ADMassFractionNuclideDepletion");
    params.set<NonlinearVariableName>("variable") = frac_var_name;
    // Set the number of groups.
    params.set<unsigned int>("num_groups") = _num_groups;

    // Set the microscopic absorption cross-sections for this isotope.
    auto & nuclide_xs = params.set<std::vector<Real>>("group_absorption");
    nuclide_xs.resize(_num_groups, 0.0);
    for (const auto & rxn : nuclide_rxns)
    {
      if (rxn._cross_sections.size() > 0u)
      {
        for (unsigned int g = 0u; g < _num_groups; ++g)
          nuclide_xs[g] += rxn._cross_sections[g];
      }
    }

    // Apply the scalar fluxes.
    auto & scalar_flux_names = params.set<std::vector<VariableName>>("group_scalar_fluxes");
    scalar_flux_names.reserve(_num_groups);
    for (unsigned int g = 0u; g < _num_groups; ++g)
      scalar_flux_names.emplace_back(_group_flux_moments[g]);

    // Apply common SUPG parameters.
    applyIsotopeParameters(params);

    if (isParamValid("block"))
    {
      params.set<std::vector<SubdomainName>>("block") =
          getParam<std::vector<SubdomainName>>("block");
    }

    _problem->addKernel(
        "ADMassFractionNuclideDepletion", "ADMassFractionNuclideDepletion" + frac_var_name, params);
    debugOutput("    - Adding kernel ADMassFractionNuclideDepletion for the variable " +
                frac_var_name + ".");
  } // ADMassFractionNuclideDepletion

  // Add ADMassFractionNuclideDiffusion.
  {
    auto params = _factory.getValidParams("ADMassFractionNuclideDiffusion");
    params.set<NonlinearVariableName>("variable") = frac_var_name;

    // Apply common SUPG parameters.
    applyIsotopeParameters(params);

    if (isParamValid("block"))
    {
      params.set<std::vector<SubdomainName>>("block") =
          getParam<std::vector<SubdomainName>>("block");
    }

    _problem->addKernel(
        "ADMassFractionNuclideDiffusion", "ADMassFractionNuclideDiffusion" + frac_var_name, params);
    debugOutput("    - Adding kernel ADMassFractionNuclideDiffusion for the variable " +
                frac_var_name + ".");
  } // ADMassFractionNuclideDiffusion

  // Add ADMassFractionNuclideTimeDerivative
  if (_exec_type == ExecutionType::Transient)
  {
    auto params = _factory.getValidParams("ADMassFractionNuclideTimeDerivative");
    params.set<NonlinearVariableName>("variable") = frac_var_name;

    // Apply common SUPG parameters.
    applyIsotopeParameters(params);

    if (isParamValid("block"))
    {
      params.set<std::vector<SubdomainName>>("block") =
          getParam<std::vector<SubdomainName>>("block");
    }

    _problem->addKernel("ADMassFractionNuclideTimeDerivative",
                        "ADMassFractionNuclideTimeDerivative" + frac_var_name,
                        params);
    debugOutput("    - Adding kernel ADMassFractionNuclideTimeDerivative for the variable " +
                frac_var_name + ".");
  } // ADMassFractionNuclideTimeDerivative
}

void
MobileDepletionSystemAction::addBCs(const std::string & nuclide_var_name)
{
  const std::string frac_var_name = nuclide_var_name + "_mass_fraction";

  // Add ADMassFractionNuclideInflowBC.
  if (_inlet_atom_fractions.size() > 0u)
  {
    // Add inlet BCs for the atom densities provided by the user.
    for (unsigned int i = 0u; i < _inlet_boundaries.size(); ++i)
    {
      auto params = _factory.getValidParams("ADMassFractionNuclideInflowBC");
      params.set<NonlinearVariableName>("variable") = frac_var_name;
      params.set<Real>("inflow_rate") = _boundary_mass_fractions[i].at(nuclide_var_name);

      applyIsotopeParameters(params);

      params.set<std::vector<BoundaryName>>("boundary").emplace_back(_inlet_boundaries[i]);

      _problem->addBoundaryCondition(
          "ADMassFractionNuclideInflowBC", _inlet_boundaries[i] + "_" + frac_var_name, params);
      debugOutput("      - Adding BC " + _inlet_boundaries[i] + "_" + frac_var_name + ".");
    }
  }
  else
  {
    // Add inlet BCs for the atom densities provided in the initial mixture.
    for (unsigned int i = 0u; i < _inlet_boundaries.size(); ++i)
    {
      auto params = _factory.getValidParams("ADMassFractionNuclideInflowBC");
      params.set<NonlinearVariableName>("variable") = frac_var_name;
      params.set<Real>("inflow_rate") = _total_nuclide_list.at(nuclide_var_name);

      applyIsotopeParameters(params);

      params.set<std::vector<BoundaryName>>("boundary").emplace_back(_inlet_boundaries[i]);

      _problem->addBoundaryCondition(
          "ADMassFractionNuclideInflowBC", _inlet_boundaries[i] + "_" + frac_var_name, params);
      debugOutput("      - Adding BC " + _inlet_boundaries[i] + "_" + frac_var_name + ".");
    }
  } // ADMassFractionNuclideInflowBC

  // Add ADMassFractionNuclideOutflowBC.
  if (_outlet_boundaries.size() > 0u)
  {
    auto params = _factory.getValidParams("ADMassFractionNuclideOutflowBC");
    params.set<NonlinearVariableName>("variable") = frac_var_name;

    applyIsotopeParameters(params);

    for (unsigned int i = 0u; i < _outlet_boundaries.size(); ++i)
      params.set<std::vector<BoundaryName>>("boundary").emplace_back(_outlet_boundaries[i]);

    _problem->addBoundaryCondition("ADMassFractionNuclideOutflowBC",
                                   "ADMassFractionNuclideOutflowBC_" + frac_var_name,
                                   params);
    debugOutput("      - Adding BC ADMassFractionNuclideOutflowBC_" + frac_var_name + ".");
  } // ADMassFractionNuclideOutflowBC
}

void
MobileDepletionSystemAction::addFVVariables(const std::string & nuclide_var_name)
{
  const std::string frac_var_name = nuclide_var_name + "_mass_fraction";

  auto params = _factory.getValidParams("INSFVScalarFieldVariable");
  params.set<std::vector<Real>>("scaling") = {getParam<Real>("scaling")};
  params.set<MooseEnum>("face_interp_method") = getParam<MooseEnum>("fv_face_interpolation");
  params.set<bool>("two_term_boundary_expansion") =
      getParam<bool>("fv_two_term_boundary_expansion");

  if (_subdomain_ids.empty())
    _problem->addVariable("INSFVScalarFieldVariable", frac_var_name, params);
  else
  {
    for (const SubdomainID & id : _subdomain_ids)
      params.set<std::vector<SubdomainName>>("block").push_back(Moose::stringify(id));

    _problem->addVariable("INSFVScalarFieldVariable", frac_var_name, params);
  }

  debugOutput("      - Adding variable " + frac_var_name + ".");
}

void
MobileDepletionSystemAction::addFVKernels(const std::string & nuclide_var_name)
{
  const std::string frac_var_name = nuclide_var_name + "_mass_fraction";

  const auto & nuclide_data = _coupled_depletion_lib->getNuclide(nuclide_var_name);
  const auto & nuclide_rxns = nuclide_data.getReactions();
  const auto & activation_parents = _coupled_depletion_lib->getActivationParents(nuclide_var_name);
  const auto & decay_parents = _coupled_depletion_lib->getDecayParents(nuclide_var_name);

  // Add ADFVMassFractionNuclideActivation.
  bool should_add_act_src = false;
  for (const auto & act_nuclide : activation_parents)
  {
    if (_total_nuclide_list.count(act_nuclide) > 0u)
    {
      should_add_act_src = true;
      break;
    }
  }
  if (should_add_act_src && _has_transport_system)
  {
    auto params = _factory.getValidParams("ADFVMassFractionNuclideActivation");
    params.set<NonlinearVariableName>("variable") = frac_var_name;
    // Set the number of groups.
    params.set<unsigned int>("num_groups") = _num_groups;

    // Set the activation cross-sections and parent mass fractions.
    auto & all_xs = params.set<std::vector<Real>>("group_activation");
    auto & fractions = params.set<std::vector<MooseFunctorName>>("isotope_mass_fractions");

    unsigned int i = 0u;
    for (const auto & act_nuclide : activation_parents)
    {
      if (_total_nuclide_list.count(act_nuclide) == 0u)
        continue;

      fractions.emplace_back(act_nuclide + "_mass_fraction");
      for (unsigned int g = 0u; g < _num_groups; ++g)
        all_xs.emplace_back(0.0);

      for (const auto & rxn : _coupled_depletion_lib->getNuclide(act_nuclide).getReactions())
      {
        if (rxn._cross_sections.size() > 0u)
        {
          for (unsigned int g = 0u; g < _num_groups; ++g)
            all_xs[i * _num_groups + g] += rxn._branching_factor * rxn._cross_sections[g];
        }
      }
      i++;
    }

    // Apply the scalar fluxes.
    auto & scalar_flux_names = params.set<std::vector<VariableName>>("group_scalar_fluxes");
    scalar_flux_names.reserve(_num_groups);
    for (unsigned int g = 0u; g < _num_groups; ++g)
      scalar_flux_names.emplace_back(_group_flux_moments[g]);

    // Apply the coupled density.
    if (_coupled_ns_fv)
      params.set<MooseFunctorName>("density") = _coupled_ns_fv->densityName();
    else
      params.set<MooseFunctorName>("density") = getParam<MooseFunctorName>("density");

    if (isParamValid("block"))
    {
      params.set<std::vector<SubdomainName>>("block") =
          getParam<std::vector<SubdomainName>>("block");
    }

    _problem->addFVKernel("ADFVMassFractionNuclideActivation",
                          "ADFVMassFractionNuclideActivation" + frac_var_name,
                          params);
    debugOutput("    - Adding kernel ADFVMassFractionNuclideActivation for the variable " +
                frac_var_name + ".");
  } // ADFVMassFractionNuclideActivation

  // Add advection kernels.
  if (_coupled_ns_fv)
  {
    // Add INSFVMassFractionScalarFieldAdvection
    {
      auto params = _factory.getValidParams("INSFVMassFractionScalarFieldAdvection");
      params.set<NonlinearVariableName>("variable") = frac_var_name;
      params.set<MooseEnum>("advected_interp_method") = getParam<MooseEnum>("fv_adv_interpolation");
      params.set<MooseEnum>("velocity_interp_method") =
          _coupled_ns_fv->getParam<MooseEnum>("velocity_interp_method");
      params.set<MooseFunctorName>("density") = _coupled_ns_fv->densityName();

      if (isParamValid("block"))
      {
        params.set<std::vector<SubdomainName>>("block") =
            getParam<std::vector<SubdomainName>>("block");
      }

      params.set<UserObjectName>("rhie_chow_user_object") = _coupled_ns_fv->rhieChowUOName();

      _problem->addFVKernel("INSFVMassFractionScalarFieldAdvection",
                            "INSFVMassFractionScalarFieldAdvection" + frac_var_name,
                            params);
      debugOutput("    - Adding kernel INSFVMassFractionScalarFieldAdvection for the variable " +
                  frac_var_name + ".");
    } // INSFVMassFractionScalarFieldAdvection
  }
  else
  {
    // Add ADFVMassFractionFunctorAdvection.
    {
      auto params = _factory.getValidParams("ADFVMassFractionFunctorAdvection");
      params.set<NonlinearVariableName>("variable") = frac_var_name;
      params.set<MooseEnum>("advected_interp_method") = getParam<MooseEnum>("fv_adv_interpolation");

      // Apply common isotope parameters.
      applyIsotopeParameters(params);

      if (isParamValid("block"))
      {
        params.set<std::vector<SubdomainName>>("block") =
            getParam<std::vector<SubdomainName>>("block");
      }

      _problem->addFVKernel("ADFVMassFractionFunctorAdvection",
                            "ADFVMassFractionFunctorAdvection" + frac_var_name,
                            params);
      debugOutput("    - Adding kernel ADFVMassFractionFunctorAdvection for the variable " +
                  frac_var_name + ".");
    } // ADFVMassFractionFunctorAdvection
  }

  // Add ADFVMassFractionNuclideDecaySink.
  if (nuclide_data.decayConst() > 0.0)
  {
    auto params = _factory.getValidParams("ADFVMassFractionNuclideDecaySink");
    params.set<NonlinearVariableName>("variable") = frac_var_name;
    // Decay constant.
    params.set<Real>("decay_const") = nuclide_data.decayConst();

    // Apply the coupled density.
    if (_coupled_ns_fv)
      params.set<MooseFunctorName>("density") = _coupled_ns_fv->densityName();
    else
      params.set<MooseFunctorName>("density") = getParam<MooseFunctorName>("density");

    if (isParamValid("block"))
    {
      params.set<std::vector<SubdomainName>>("block") =
          getParam<std::vector<SubdomainName>>("block");
    }

    _problem->addFVKernel("ADFVMassFractionNuclideDecaySink",
                          "ADFVMassFractionNuclideDecaySink" + frac_var_name,
                          params);
    debugOutput("    - Adding kernel ADFVMassFractionNuclideDecaySink for the variable " +
                frac_var_name + ".");
  } // ADFVMassFractionNuclideDecaySink

  // Add ADFVMassFractionNuclideDecaySource.
  bool should_add_decay_src = false;
  for (const auto & dec_nuclide : decay_parents)
  {
    if (_total_nuclide_list.count(dec_nuclide) > 0u)
    {
      should_add_decay_src = true;
      break;
    }
  }
  if (should_add_decay_src)
  {
    auto params = _factory.getValidParams("ADFVMassFractionNuclideDecaySource");
    params.set<NonlinearVariableName>("variable") = frac_var_name;

    // Apply decay constants and branching factors.
    auto & decay_data = params.set<std::vector<Real>>("decay_constants");

    auto & fractions = params.set<std::vector<MooseFunctorName>>("isotope_fractions");

    unsigned int i = 0u;
    for (const auto & dec_nuclide : decay_parents)
    {
      if (_total_nuclide_list.count(dec_nuclide) == 0u)
        continue;

      fractions.emplace_back(dec_nuclide + "_mass_fraction");
      decay_data.emplace_back(0.0);

      for (const auto & decay : _coupled_depletion_lib->getNuclide(dec_nuclide).getDecays())
      {
        if (decay._target == nuclide_var_name)
          decay_data[i] += decay._branching_factor;
      }
      decay_data[i] *= _coupled_depletion_lib->getNuclide(dec_nuclide).decayConst();
      i++;
    }

    // Apply the coupled density.
    if (_coupled_ns_fv)
      params.set<MooseFunctorName>("density") = _coupled_ns_fv->densityName();
    else
      params.set<MooseFunctorName>("density") = getParam<MooseFunctorName>("density");

    if (isParamValid("block"))
    {
      params.set<std::vector<SubdomainName>>("block") =
          getParam<std::vector<SubdomainName>>("block");
    }

    _problem->addFVKernel("ADFVMassFractionNuclideDecaySource",
                          "ADFVMassFractionNuclideDecaySource" + frac_var_name,
                          params);
    debugOutput("    - Adding kernel ADFVMassFractionNuclideDecaySource for the variable " +
                frac_var_name + ".");
  } // ADFVMassFractionNuclideDecaySource

  // Add ADFVMassFractionNuclideDepletion.
  if (nuclide_rxns.size() > 0u && _has_transport_system)
  {
    auto params = _factory.getValidParams("ADFVMassFractionNuclideDepletion");
    params.set<NonlinearVariableName>("variable") = frac_var_name;
    // Set the number of groups.
    params.set<unsigned int>("num_groups") = _num_groups;

    // Set the microscopic absorption cross-sections for this isotope.
    auto & nuclide_xs = params.set<std::vector<Real>>("group_absorption");
    nuclide_xs.resize(_num_groups, 0.0);
    for (const auto & rxn : nuclide_rxns)
    {
      if (rxn._cross_sections.size() > 0u)
      {
        for (unsigned int g = 0u; g < _num_groups; ++g)
          nuclide_xs[g] += rxn._cross_sections[g];
      }
    }

    // Apply the scalar fluxes.
    auto & scalar_flux_names = params.set<std::vector<VariableName>>("group_scalar_fluxes");
    scalar_flux_names.reserve(_num_groups);
    for (unsigned int g = 0u; g < _num_groups; ++g)
      scalar_flux_names.emplace_back(_group_flux_moments[g]);

    // Apply the coupled density.
    if (_coupled_ns_fv)
      params.set<MooseFunctorName>("density") = _coupled_ns_fv->densityName();
    else
      params.set<MooseFunctorName>("density") = getParam<MooseFunctorName>("density");

    if (isParamValid("block"))
    {
      params.set<std::vector<SubdomainName>>("block") =
          getParam<std::vector<SubdomainName>>("block");
    }

    _problem->addFVKernel("ADFVMassFractionNuclideDepletion",
                          "ADFVMassFractionNuclideDepletion" + frac_var_name,
                          params);
    debugOutput("    - Adding kernel ADFVMassFractionNuclideDepletion for the variable " +
                frac_var_name + ".");
  } // ADFVMassFractionNuclideDepletion

  // Add ADFVMassFractionNuclideDiffusion.
  {
    auto params = _factory.getValidParams("ADFVMassFractionNuclideDiffusion");
    params.set<NonlinearVariableName>("variable") = frac_var_name;

    // Apply the coupled density.
    if (_coupled_ns_fv)
      params.set<MooseFunctorName>("density") = _coupled_ns_fv->densityName();
    else
      params.set<MooseFunctorName>("density") = getParam<MooseFunctorName>("density");

    if (isParamValid("block"))
    {
      params.set<std::vector<SubdomainName>>("block") =
          getParam<std::vector<SubdomainName>>("block");
    }

    _problem->addFVKernel("ADFVMassFractionNuclideDiffusion",
                          "ADFVMassFractionNuclideDiffusion" + frac_var_name,
                          params);
    debugOutput("    - Adding kernel ADFVMassFractionNuclideDiffusion for the variable " +
                frac_var_name + ".");
  } // ADFVMassFractionNuclideDiffusion

  // Add ADFVMassFractionTimeDerivative
  if (_exec_type == ExecutionType::Transient)
  {
    auto params = _factory.getValidParams("ADFVMassFractionTimeDerivative");
    params.set<NonlinearVariableName>("variable") = frac_var_name;

    // Apply the coupled density.
    if (_coupled_ns_fv)
      params.set<MooseFunctorName>("density") = _coupled_ns_fv->densityName();
    else
      params.set<MooseFunctorName>("density") = getParam<MooseFunctorName>("density");

    if (isParamValid("block"))
    {
      params.set<std::vector<SubdomainName>>("block") =
          getParam<std::vector<SubdomainName>>("block");
    }

    _problem->addFVKernel(
        "ADFVMassFractionTimeDerivative", "ADFVMassFractionTimeDerivative" + frac_var_name, params);
    debugOutput("    - Adding kernel ADFVMassFractionTimeDerivative for the variable " +
                frac_var_name + ".");
  } // ADFVMassFractionTimeDerivative
}

void
MobileDepletionSystemAction::addFVBCs(const std::string & nuclide_var_name)
{
  const std::string frac_var_name = nuclide_var_name + "_mass_fraction";

  // Add ADFVMassFractionNuclideInflowBC
  if (_inlet_atom_fractions.size() > 0u)
  {
    // Add inlet BCs for the atom densities provided by the user.
    for (unsigned int i = 0u; i < _inlet_boundaries.size(); ++i)
    {
      auto params = _factory.getValidParams("ADFVMassFractionNuclideInflowBC");
      params.set<NonlinearVariableName>("variable") = frac_var_name;
      params.set<Real>("inflow_rate") = _boundary_mass_fractions[i].at(nuclide_var_name);

      applyIsotopeParameters(params);

      params.set<std::vector<BoundaryName>>("boundary").emplace_back(_inlet_boundaries[i]);

      _problem->addFVBC(
          "ADFVMassFractionNuclideInflowBC", _inlet_boundaries[i] + "_" + frac_var_name, params);
      debugOutput("      - Adding FVBC " + _inlet_boundaries[i] + "_" + frac_var_name + ".");
    }
  }
  else
  {
    // Add inlet BCs for the atom densities provided in the initial mixture.
    for (unsigned int i = 0u; i < _inlet_boundaries.size(); ++i)
    {
      auto params = _factory.getValidParams("ADFVMassFractionNuclideInflowBC");
      params.set<NonlinearVariableName>("variable") = frac_var_name;
      params.set<Real>("inflow_rate") = _total_nuclide_list.at(nuclide_var_name);

      applyIsotopeParameters(params);

      params.set<std::vector<BoundaryName>>("boundary").emplace_back(_inlet_boundaries[i]);

      _problem->addFVBC(
          "ADFVMassFractionNuclideInflowBC", _inlet_boundaries[i] + "_" + frac_var_name, params);
      debugOutput("      - Adding FVBC " + _inlet_boundaries[i] + "_" + frac_var_name + ".");
    }
  } // ADFVMassFractionNuclideInflowBC
}

void
MobileDepletionSystemAction::addAuxVariables(const std::string & nuclide_var_name)
{
  if (_scheme == NuclideScheme::FV)
  {
    auto params = _factory.getValidParams("MooseVariableFVReal");
    params.set<bool>("two_term_boundary_expansion") =
        getParam<bool>("fv_two_term_boundary_expansion");

    if (isParamValid("block"))
    {
      params.set<std::vector<SubdomainName>>("block") =
          getParam<std::vector<SubdomainName>>("block");
    }

    _problem->addAuxVariable("MooseVariableFVReal", nuclide_var_name, params);
    debugOutput("      - Adding auxvariable " + nuclide_var_name + ".");
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

    _problem->addAuxVariable(type, nuclide_var_name, params);
    debugOutput("      - Adding auxvariable " + nuclide_var_name + ".");
  }
}

void
MobileDepletionSystemAction::addAuxKernels(const std::string & nuclide_var_name)
{
  const std::string frac_var_name = nuclide_var_name + "_mass_fraction";

  // Add MassFractionConcentration
  {
    auto params = _factory.getValidParams("MassFractionConcentration");
    params.set<AuxVariableName>("variable") = nuclide_var_name;
    params.set<bool>("compute_number_density") = getParam<bool>("compute_number_density");

    if (_scheme == NuclideScheme::SUPGFE)
    {
      params.set<std::vector<VariableName>>("nuclide_mass_fraction_var")
          .emplace_back(frac_var_name);
      params.set<bool>("is_fe") = true;
    }
    else
    {
      params.set<MooseFunctorName>("nuclide_mass_fraction_fun") = frac_var_name;
      params.set<bool>("is_fe") = false;
    }

    // Set the nuclide name.
    params.set<std::string>("nuclide_name") = nuclide_var_name;

    // Apply the coupled density.
    if (_coupled_ns_fv)
      params.set<MooseFunctorName>("density") = _coupled_ns_fv->densityName();
    else
      params.set<MooseFunctorName>("density") = getParam<MooseFunctorName>("density");

    if (isParamValid("block"))
    {
      params.set<std::vector<SubdomainName>>("block") =
          getParam<std::vector<SubdomainName>>("block");
    }

    params.set<ExecFlagEnum>("execute_on") = {
        EXEC_INITIAL, EXEC_TIMESTEP_BEGIN, EXEC_TIMESTEP_END, EXEC_LINEAR};

    _problem->addAuxKernel(
        "MassFractionConcentration", "MassFractionConcentration_" + nuclide_var_name, params);
    debugOutput("      - Adding auxkernel MassFractionConcentration for the variable " +
                nuclide_var_name + ".");
  } // MassFractionConcentration
}

// Fetch properties from the coupled TransportSystem.
void
MobileDepletionSystemAction::fetchTransportProperties()
{
  if (_transport_system == "")
  {
    // Doesn't add any kernels related to neutron activation / depletion.
    _has_transport_system = false;
    return;
  }

  // Fetches the required parameters from the coupled TransportSystem and prepares variable names.
  const auto & transport_action = _awh.getAction<TransportAction>(_transport_system);
  _num_groups = transport_action.getParam<unsigned int>("num_groups");
  const std::string & flux_moment_name =
      transport_action.getParam<std::string>("flux_moment_names");
  for (unsigned int g = 0; g < _num_groups; ++g)
  {
    // Set up variable names for the group flux moments.
    _group_flux_moments.emplace_back(flux_moment_name + "_" + Moose::stringify(g + 1u) + "_" +
                                     Moose::stringify(0) + "_" + Moose::stringify(0));
  }
}

// Modify outputs to remove the mass fraction primal variables.
void
MobileDepletionSystemAction::modifyOutputs()
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
      for (const auto & [nuclide, fraction] : _total_nuclide_list)
        output_params.set<std::vector<VariableName>>("hide").emplace_back(nuclide +
                                                                          "_mass_fraction");
    }
  }
}

// Auxvariables and auxkernels for radiation transport source terms.
void
MobileDepletionSystemAction::addRadiationAuxVariables()
{
  if (getParam<bool>("add_photon_sources"))
  {
    for (unsigned int g = 0u; g < _photon_group_boundaries.size() - 1u; ++g)
    {
      if (_scheme == NuclideScheme::FV)
      {
        auto params = _factory.getValidParams("MooseVariableFVReal");
        params.set<bool>("two_term_boundary_expansion") =
            getParam<bool>("fv_two_term_boundary_expansion");

        if (isParamValid("block"))
        {
          params.set<std::vector<SubdomainName>>("block") =
              getParam<std::vector<SubdomainName>>("block");
        }

        _problem->addAuxVariable(
            "MooseVariableFVReal", _photon_source_prefix + "_g_" + Moose::stringify(g), params);
        debugOutput("      - Adding auxvariable " + _photon_source_prefix + "_g_" +
                    Moose::stringify(g) + ".");
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

        _problem->addAuxVariable(type, _photon_source_prefix + "_g_" + Moose::stringify(g), params);
        debugOutput("      - Adding auxvariable " + _photon_source_prefix + "_g_" +
                    Moose::stringify(g) + ".");
      }
    }
  }

  if (getParam<bool>("add_neutron_sources"))
  {
    for (unsigned int g = 0u; g < _neutron_group_boundaries.size() - 1u; ++g)
    {
      if (_scheme == NuclideScheme::FV)
      {
        auto params = _factory.getValidParams("MooseVariableFVReal");
        params.set<bool>("two_term_boundary_expansion") =
            getParam<bool>("fv_two_term_boundary_expansion");

        if (isParamValid("block"))
        {
          params.set<std::vector<SubdomainName>>("block") =
              getParam<std::vector<SubdomainName>>("block");
        }

        _problem->addAuxVariable(
            "MooseVariableFVReal", _neutron_source_prefix + "_g_" + Moose::stringify(g), params);
        debugOutput("      - Adding auxvariable " + _neutron_source_prefix + "_g_" +
                    Moose::stringify(g) + ".");
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

        _problem->addAuxVariable(
            type, _neutron_source_prefix + "_g_" + Moose::stringify(g), params);
        debugOutput("      - Adding auxvariable " + _neutron_source_prefix + "_g_" +
                    Moose::stringify(g) + ".");
      }
    }
  }
}

void
MobileDepletionSystemAction::addRadiationAuxKernels()
{
  if (getParam<bool>("add_photon_sources"))
  {
    for (unsigned int g = 0u; g < _photon_group_boundaries.size() - 1u; ++g)
    {
      auto params = _factory.getValidParams("ParticleDecaySource");
      params.set<AuxVariableName>("variable") = _photon_source_prefix + "_g_" + Moose::stringify(g);
      params.set<MooseEnum>("particle_type") = MooseEnum("neutron photon", "photon");
      params.set<std::vector<Real>>("group_boundaries") = _photon_group_boundaries;
      params.set<unsigned int>("group_index") = g;
      params.set<bool>("is_number_density") = getParam<bool>("compute_number_density");

      if (_scheme == NuclideScheme::FV)
      {
        params.set<bool>("is_fe") = false;
        for (const auto & [nuclide, weight_fraction] : _total_nuclide_list)
          params.set<std::vector<MooseFunctorName>>("nuclide_funs").emplace_back(nuclide);
      }
      else
      {
        params.set<bool>("is_fe") = true;
        for (const auto & [nuclide, weight_fraction] : _total_nuclide_list)
          params.set<std::vector<VariableName>>("nuclide_vars").emplace_back(nuclide);
      }

      params.set<UserObjectName>("data_lib_name") =
          _coupled_depletion_lib->getParam<std::string>("depletion_uo_name");

      if (isParamValid("block"))
      {
        params.set<std::vector<SubdomainName>>("block") =
            getParam<std::vector<SubdomainName>>("block");
      }

      _problem->addAuxKernel(
          "ParticleDecaySource", "ParticleDecaySource_Photon_" + Moose::stringify(g), params);
      debugOutput("      - Adding auxkernel ParticleDecaySource for the variable " +
                  _photon_source_prefix + "_g_" + Moose::stringify(g) + ".");
    }
  }

  if (getParam<bool>("add_neutron_sources"))
  {
    for (unsigned int g = 0u; g < _neutron_group_boundaries.size() - 1u; ++g)
    {
      auto params = _factory.getValidParams("ParticleDecaySource");
      params.set<AuxVariableName>("variable") =
          _neutron_source_prefix + "_g_" + Moose::stringify(g);
      params.set<MooseEnum>("particle_type") = MooseEnum("neutron photon", "neutron");
      params.set<std::vector<Real>>("group_boundaries") = _neutron_group_boundaries;
      params.set<unsigned int>("group_index") = g;
      params.set<bool>("is_number_density") = getParam<bool>("compute_number_density");

      if (_scheme == NuclideScheme::FV)
      {
        params.set<bool>("is_fe") = false;
        for (const auto & [nuclide, weight_fraction] : _total_nuclide_list)
          params.set<std::vector<MooseFunctorName>>("nuclide_funs").emplace_back(nuclide);
      }
      else
      {
        params.set<bool>("is_fe") = true;
        for (const auto & [nuclide, weight_fraction] : _total_nuclide_list)
          params.set<std::vector<VariableName>>("nuclide_vars").emplace_back(nuclide);
      }

      params.set<UserObjectName>("data_lib_name") =
          _coupled_depletion_lib->getParam<std::string>("depletion_uo_name");

      if (isParamValid("block"))
      {
        params.set<std::vector<SubdomainName>>("block") =
            getParam<std::vector<SubdomainName>>("block");
      }

      _problem->addAuxKernel(
          "ParticleDecaySource", "ParticleDecaySource_Neutron_" + Moose::stringify(g), params);
      debugOutput("      - Adding auxkernel ParticleDecaySource for the variable " +
                  _neutron_source_prefix + "_g_" + Moose::stringify(g) + ".");
    }
  }
}

void
MobileDepletionSystemAction::act()
{
  if (_first_action)
  {
    debugOutput("Initializing a mass transport system: ", "Initializing a mass transport system: ");

    initializeBase();

    // Fetch the required coupled depletion library.
    {
      const auto depletion_lib_actions = _awh.getActions<DepletionLibraryAction>();
      if (depletion_lib_actions.size() != 1u)
        mooseError("The input file must have a single Depletion Library. There are currently " +
                   Moose::stringify(depletion_lib_actions.size() + " DepletionLibraryActions."));
      _coupled_depletion_lib = depletion_lib_actions[0u];
    }

    // Fetch the coupled Navier-Stokes finite volume physics.
    if (_using_moose_ns)
    {
      const auto ns_fv_physics = _awh.getActions<WCNSFVFlowPhysics>();
      if (ns_fv_physics.size() != 1u)
        mooseError("The input file must have a single Navier-Stokes system. There are currently " +
                   Moose::stringify(ns_fv_physics.size() + " WCNSFVFlowPhysics."));
      _coupled_ns_fv = ns_fv_physics[0u];
    }

    // Fetch the properties from the transport system.
    fetchTransportProperties();
    // Add all elements.
    addElementsAndNuclides();
    // Fetch depletion properties and build the depletion chain.
    buildDepletionSystem();

    _mesh_dims = _problem->mesh().dimension();
    _first_action = false;
  }

  // Modify outputs to only output the resulting density auxvariables.
  if (_current_task == "add_variable" && !getParam<bool>("output_nuclide_mass_fractions"))
  {
    debugOutput("    - Modifying outputs...");
    modifyOutputs();
  }

  for (const auto & [nuclide, weight_fraction] : _total_nuclide_list)
  {
    if (_current_task == "add_variable" && _scheme == NuclideScheme::SUPGFE)
    {
      debugOutput("  - Adding variables...");

      addVariable(nuclide + "_mass_fraction");
    }

    if (_current_task == "add_variable" && _scheme == NuclideScheme::FV)
    {
      debugOutput("  - Adding variables...");
      addFVVariables(nuclide);
    }

    if (_current_task == "add_ic")
    {
      debugOutput("  - Adding initial conditions...");
      addICs(nuclide);
    }

    if (_current_task == "add_kernel" && _scheme == NuclideScheme::SUPGFE)
    {
      debugOutput("  - Adding kernels...");
      addKernels(nuclide);
    }

    if (_current_task == "add_bc" && _scheme == NuclideScheme::SUPGFE)
    {
      debugOutput("  - Adding BCs...");
      addBCs(nuclide);
    }

    if (_current_task == "add_fv_kernel" && _scheme == NuclideScheme::FV)
    {
      debugOutput("  - Adding FV kernels...");
      addFVKernels(nuclide);
    }

    if (_current_task == "add_fv_bc" && _scheme == NuclideScheme::FV)
    {
      debugOutput("  - Adding FV BCs...");
      addFVBCs(nuclide);
    }

    if (_current_task == "add_material")
    {
      debugOutput("  - Adding materials...");
      addMaterials(nuclide);
    }

    if (_current_task == "add_aux_variable")
    {
      debugOutput("  - Adding auxvariables...");
      addAuxVariables(nuclide);
    }

    if (_current_task == "add_aux_kernel")
    {
      debugOutput("  - Adding auxkernels...");
      addAuxKernels(nuclide);
    }
  }

  // Add variables and auxkernels for radionuclide particle sources.
  if (getParam<bool>("add_photon_sources") || getParam<bool>("add_neutron_sources"))
  {
    if (_current_task == "add_aux_variable")
    {
      debugOutput("  - Adding radiation transport auxvariables...");
      addRadiationAuxVariables();
    }

    if (_current_task == "add_aux_kernel")
    {
      debugOutput("  - Adding radiation transport auxkernels...");
      addRadiationAuxKernels();
    }
  }
}
