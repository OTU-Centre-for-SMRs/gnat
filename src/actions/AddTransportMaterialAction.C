#include "AddTransportMaterialAction.h"

#include "ActionWarehouse.h"
#include "FEProblemBase.h"

#include "TransportAction.h"

registerMooseAction("GnatApp", AddTransportMaterialAction, "add_material");

InputParameters
AddTransportMaterialAction::validParams()
{
  auto params = MooseObjectAction::validParams();
  params.addClassDescription(
      "Adds a material using provided parameters and information from a transport system.");
  params.addParam<std::string>(
      "transport_system",
      "",
      "Name of the transport system which will consume the provided material properties. If one is "
      "not provided the first transport system will be used.");

  return params;
}

AddTransportMaterialAction::AddTransportMaterialAction(const InputParameters & parameters)
  : MooseObjectAction(parameters),
    _parent_transport_system(getParam<std::string>("transport_system")),
    _num_groups(0u),
    _particle(MooseEnum("neutron photon", "neutron")),
    _scheme(MooseEnum("saaf_cfem diffusion_cfem flux_moment_transfer")),
    _disable_fission(true),
    _is_init(false)
{
}

void
AddTransportMaterialAction::act()
{
  if (!_is_init)
  {
    // Get the associated transport action.
    if (_parent_transport_system == "")
    {
      // A transport system was not provided. Fetch the first transport system from the action
      // warehouse.
      auto transport_actions = _awh.getActions<TransportAction>();
      if (transport_actions.size() > 0u)
      {
        _num_groups = transport_actions[0u]->getParam<unsigned int>("num_groups");
        _particle = transport_actions[0u]->getParam<MooseEnum>("particle_type");
        _scheme = transport_actions[0u]->getParam<MooseEnum>("scheme");
        _disable_fission = transport_actions[0u]->getParam<bool>("debug_disable_fission");

        _parent_transport_system = transport_actions[0u]->name();
      }
      else if (transport_actions.size() == 0u)
        mooseError("No transport systems have been declared in the input deck.");
      else
        mooseError(
            "Multiple transport systems have been declared. Please select one in the input deck.");
    }
    else
    {
      const auto & transport_action = _awh.getAction<TransportAction>(_parent_transport_system);
      _num_groups = transport_action.getParam<unsigned int>("num_groups");
      _particle = transport_action.getParam<MooseEnum>("particle_type");
      _scheme = transport_action.getParam<MooseEnum>("scheme");
      _disable_fission = transport_action.getParam<bool>("debug_disable_fission");
    }

    _moose_object_pars.set<unsigned int>("num_groups") = _num_groups;
    _moose_object_pars.set<MooseEnum>("particle_type") = _particle;
    _moose_object_pars.set<bool>("is_saaf") = _scheme == "saaf_cfem";
    _moose_object_pars.set<bool>("is_diffusion") = _scheme == "diffusion_cfem";
    _moose_object_pars.set<bool>("has_fission") = _particle == "neutron" && !_disable_fission;
    _moose_object_pars.set<std::string>("transport_system") = _parent_transport_system;

    _is_init = true;
  }

  _problem->addMaterial(_type, _name, _moose_object_pars);
}
