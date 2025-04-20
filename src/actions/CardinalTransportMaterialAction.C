#include "CardinalTransportMaterialAction.h"

#include "TransportAction.h"
#include "UncollidedFluxAction.h"

#include "AuxiliarySystem.h"

registerMooseAction("GnatApp", CardinalTransportMaterialAction, "add_aux_variable");
registerMooseAction("GnatApp", CardinalTransportMaterialAction, "add_material");

registerMooseAction("GnatApp", CardinalTransportMaterialAction, "check_copy_nodal_vars");
registerMooseAction("GnatApp", CardinalTransportMaterialAction, "copy_nodal_vars");

// TODO: add transfers to pull from a Cardinal sub-app.

InputParameters
CardinalTransportMaterialAction::validParams()
{
  auto params = GnatBaseAction::validParams();
  params.addClassDescription("An action which automates the addition of Cardinal-generated cross "
                             "sections to a Gnat problem.");
  params.addRequiredParam<MooseEnum>(
      "xs_source",
      MooseEnum("mesh subapp", "mesh"),
      "The source of the cross sections. These can either be the mesh (through an exodus restart) "
      "or a Cardinal sub-application.");

  params.addParam<std::string>(
      "transport_system",
      "",
      "The name of the transport system the created material will be supplying properties to. If "
      "not provided, this action will choose the first transport system added to the problem.");

  params.addRequiredParam<unsigned int>("scatter_anisotropy",
                                        "The maximum degree of scattering anisotropy.");

  params.addParam<bool>(
      "add_fission_heating",
      false,
      "Whether or not per-group fission heating (kappa-fission) values should be generated.");

  return params;
}

CardinalTransportMaterialAction::CardinalTransportMaterialAction(const InputParameters & parameters)
  : GnatBaseAction(parameters),
    _xs_source(getParam<MooseEnum>("xs_source")),
    _parent_transport_system(getParam<std::string>("transport_system")),
    _particle(MooseEnum("neutron photon", "neutron")),
    _scheme(MooseEnum("saaf_cfem diffusion_cfem flux_moment_transfer")),
    _anisotropy(getParam<unsigned int>("scatter_anisotropy")),
    _is_init(false),
    _add_kappa_fission(getParam<bool>("add_fission_heating"))
{
}

void
CardinalTransportMaterialAction::act()
{
  if (!_is_init && _problem)
  {
    // Get the associated transport action.
    auto transport_actions = _awh.getActions<TransportAction>();
    auto uncollided_flux_actions = _awh.getActions<UncollidedFluxAction>();
    bool can_setup = true;
    if (_parent_transport_system == "")
    {
      // A transport system was not provided. Fetch the first transport system from the action
      // warehouse. Prioritizing transport systems first.
      if (transport_actions.size() == 1u)
      {
        _num_groups = transport_actions[0u]->getParam<unsigned int>("num_groups");
        _particle = transport_actions[0u]->getParam<MooseEnum>("particle_type");
        _scheme = transport_actions[0u]->getParam<MooseEnum>("scheme");
        _disable_fission = transport_actions[0u]->getParam<bool>("debug_disable_fission");

        _parent_transport_system = transport_actions[0u]->name();

        can_setup = false;
      }

      if (can_setup && uncollided_flux_actions.size() == 1u)
      {
        _num_groups = uncollided_flux_actions[0u]->getParam<unsigned int>("num_groups");
        _particle = MooseEnum("neutron photon", "photon");
        _scheme = MooseEnum("saaf_cfem diffusion_cfem flux_moment_transfer", "saaf_cfem");
        _disable_fission = true;

        _parent_transport_system = uncollided_flux_actions[0u]->name();

        can_setup = false;
      }

      if (transport_actions.size() == 0u && uncollided_flux_actions.size() == 0u)
        mooseError(
            "No transport systems / uncollided flux systems have been declared in the input deck.");

      if (transport_actions.size() > 1u || uncollided_flux_actions.size() > 1u)
        mooseError("Multiple transport systems / uncollided flux systems have been declared. "
                   "Please select one in the input deck.");
    }
    else
    {
      if (transport_actions.size() >= 1u)
      {
        const auto & transport_action = _awh.getAction<TransportAction>(_parent_transport_system);
        _num_groups = transport_action.getParam<unsigned int>("num_groups");
        _particle = transport_action.getParam<MooseEnum>("particle_type");
        _scheme = transport_action.getParam<MooseEnum>("scheme");
        _disable_fission = transport_action.getParam<bool>("debug_disable_fission");

        can_setup = false;
      }

      if (can_setup && uncollided_flux_actions.size() >= 1u)
      {
        const auto & uncollided_flux_action =
            _awh.getAction<UncollidedFluxAction>(_parent_transport_system);
        _num_groups = uncollided_flux_action.getParam<unsigned int>("num_groups");
        _particle = MooseEnum("neutron photon", "photon");
        _scheme = MooseEnum("saaf_cfem diffusion_cfem flux_moment_transfer", "saaf_cfem");
        _disable_fission = true;

        _parent_transport_system = uncollided_flux_action.name();

        can_setup = false;
      }
    }

    for (unsigned int g = 0; g < _num_groups; ++g)
    {
      _total_var_names.emplace_back("total_xs_g" + Moose::stringify(g + 1u));

      for (unsigned int gp = 0u; gp < _num_groups; ++gp)
        for (unsigned int l = 0u; l <= _anisotropy; ++l)
          _scattering_var_names.emplace_back("scatter_xs_g" + Moose::stringify(g + 1u) + "_gp" +
                                             Moose::stringify(gp + 1u) + "_l" +
                                             Moose::stringify(l));

      if (_particle == "neutron" && !_disable_fission)
      {
        _nu_fission_var_names.emplace_back("nu_fission_xs_g" + Moose::stringify(g + 1u));
        _chi_var_names.emplace_back("chi_g" + Moose::stringify(g + 1u));

        if (_add_kappa_fission)
          _kf_var_names.emplace_back("kappa_fission_g" + Moose::stringify(g + 1u));
      }

      if (_problem->isTransient())
        _inv_v_var_names.emplace_back("inv_v_g" + Moose::stringify(g + 1u));

      if (_scheme == "diffusion_cfem")
      {
        _diff_var_names.emplace_back("diff_g" + Moose::stringify(g + 1u));
        _abs_var_names.emplace_back("abs_xs_g" + Moose::stringify(g + 1u));
      }
    }

    // Initialize the base action.
    initializeBase();

    _is_init = true;
  }

  if (_current_task == "add_aux_variable")
    addAuxVariables();

  if (_current_task == "add_material")
    addMaterials();

  if (_current_task == "check_copy_nodal_vars" && _xs_source == "mesh")
    _app.setExodusFileRestart(true);

  if (_current_task == "copy_nodal_vars" && _xs_source == "mesh")
    copyOnRestart();
}

void
CardinalTransportMaterialAction::addAuxVariables()
{
  for (const auto & n : _total_var_names)
    addMaterialAuxVar(n);

  for (const auto & n : _scattering_var_names)
    addMaterialAuxVar(n);

  for (const auto & n : _nu_fission_var_names)
    addMaterialAuxVar(n);

  for (const auto & n : _chi_var_names)
    addMaterialAuxVar(n);

  for (const auto & n : _kf_var_names)
    addMaterialAuxVar(n);

  for (const auto & n : _inv_v_var_names)
    addMaterialAuxVar(n);

  for (const auto & n : _diff_var_names)
    addMaterialAuxVar(n);

  for (const auto & n : _abs_var_names)
    addMaterialAuxVar(n);
}

void
CardinalTransportMaterialAction::addMaterialAuxVar(const std::string & name)
{
  auto params = _factory.getValidParams("MooseVariable");
  params.set<MooseEnum>("family") = "MONOMIAL";
  params.set<MooseEnum>("order") = "CONSTANT";

  for (const SubdomainID & id : _subdomain_ids)
    params.set<std::vector<SubdomainName>>("block").push_back(Moose::stringify(id));

  _problem->addAuxVariable("MooseVariable", name, params);
}

void
CardinalTransportMaterialAction::addMaterials()
{
  auto params = _factory.getValidParams("PropsFromVarTransportMaterial");

  params.set<unsigned int>("num_groups") = _num_groups;
  params.set<MooseEnum>("particle_type") = _particle;
  params.set<bool>("is_saaf") = _scheme == "saaf_cfem";
  params.set<bool>("is_diffusion") = _scheme == "diffusion_cfem";
  params.set<bool>("has_fission") = _particle == "neutron" && !_disable_fission;
  params.set<std::string>("transport_system") = _parent_transport_system;
  params.set<bool>("add_heating") = _add_kappa_fission;

  for (const auto & n : _total_var_names)
    params.set<std::vector<VariableName>>("total_xs").emplace_back(n);

  for (const auto & n : _scattering_var_names)
    params.set<std::vector<VariableName>>("scatter_xs").emplace_back(n);

  for (const auto & n : _nu_fission_var_names)
    params.set<std::vector<VariableName>>("nu_fission_xs").emplace_back(n);

  for (const auto & n : _chi_var_names)
    params.set<std::vector<VariableName>>("chi").emplace_back(n);

  for (const auto & n : _kf_var_names)
    params.set<std::vector<VariableName>>("kappa_fission").emplace_back(n);

  for (const auto & n : _inv_v_var_names)
    params.set<std::vector<VariableName>>("inv_vel").emplace_back(n);

  for (const auto & n : _diff_var_names)
    params.set<std::vector<VariableName>>("diffusion").emplace_back(n);

  for (const auto & n : _abs_var_names)
    params.set<std::vector<VariableName>>("absorption_xs").emplace_back(n);

  for (const SubdomainID & id : _subdomain_ids)
    params.set<std::vector<SubdomainName>>("block").push_back(Moose::stringify(id));

  _problem->addMaterial(
      "PropsFromVarTransportMaterial", "CardinalXS_" + _parent_transport_system, params);
}

void
CardinalTransportMaterialAction::copyOnRestart()
{
  auto & aux_system = _problem->getAuxiliarySystem();

  for (const auto & n : _total_var_names)
    aux_system.addVariableToCopy(n, n, "LATEST");

  for (const auto & n : _scattering_var_names)
    aux_system.addVariableToCopy(n, n, "LATEST");

  for (const auto & n : _nu_fission_var_names)
    aux_system.addVariableToCopy(n, n, "LATEST");

  for (const auto & n : _chi_var_names)
    aux_system.addVariableToCopy(n, n, "LATEST");

  for (const auto & n : _kf_var_names)
    aux_system.addVariableToCopy(n, n, "LATEST");

  for (const auto & n : _inv_v_var_names)
    aux_system.addVariableToCopy(n, n, "LATEST");

  for (const auto & n : _diff_var_names)
    aux_system.addVariableToCopy(n, n, "LATEST");

  for (const auto & n : _abs_var_names)
    aux_system.addVariableToCopy(n, n, "LATEST");
}
