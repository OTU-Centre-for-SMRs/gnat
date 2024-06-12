#include "DepletionLibraryAction.h"

#include <queue>
#include <algorithm>

#include <filesystem>
#include "pugixml.h"

#include "Factory.h"
#include "FEProblemBase.h"

registerMooseAction("GnatApp", DepletionLibraryAction, "meta_action");
registerMooseAction("GnatApp", DepletionLibraryAction, "add_user_object");

InputParameters
DepletionLibraryAction::validParams()
{
  auto params = Action::validParams();
  params.addClassDescription("An action which provides depletion data to other MOOSE objects.");
  params.addRequiredParam<std::string>("depletion_file", "The name of the depletion file.");
  params.addParam<std::string>(
      "cross_section_file", "", "The name of the microscopic depletion cross-section file.");
  params.addParam<MooseEnum>("depletion_file_source",
                             MooseEnum("openmc", "openmc"),
                             "The depletion file generator. More depletion file sources will be "
                             "supported at a later date.");
  params.addParam<bool>(
      "show_warnings",
      false,
      "Whether or not the depletion library should output warnings to the console.");
  params.addParam<bool>(
      "add_data_userobject",
      true,
      "Whether or not the depletion data provider user object should be added to the simulation.");
  params.addParam<std::string>(
      "depletion_uo_name", "DepletionDataProviderUO", "The name of the depletion user object.");

  return params;
}

DepletionLibraryAction::DepletionLibraryAction(const InputParameters & parameters)
  : Action(parameters),
    _file_source(getParam<MooseEnum>("depletion_file_source").getEnum<DepletionFileSource>()),
    _depletion_file_name((std::filesystem::path(_app.getLastInputFileName()).parent_path() /
                          std::filesystem::path(getParam<std::string>("depletion_file")))
                             .string()),
    _mico_xs_file_name(getParam<std::string>("cross_section_file") != ""
                           ? (std::filesystem::path(_app.getLastInputFileName()).parent_path() /
                              std::filesystem::path(getParam<std::string>("cross_section_file")))
                                 .string()
                           : ""),
    _warnings(getParam<bool>("show_warnings")),
    _xs_units(XSUnits::Barns),
    _convert_to_eV(1.0),
    _num_groups(1u)
{
  _console << COLOR_YELLOW
           << "Warning: GNAT does not track secondary particles produced during radioactive decay "
              "or during neutron reactions (i.e, alpha particles in alpha decays or (n,alpha) "
              "reactions).\n"
           << COLOR_DEFAULT;
  switch (_file_source)
  {
    case DepletionFileSource::OpenMC:
      _console << "Parsing depletion library..." << std::endl;
      loadOpenMCDepletionXML();
      _console << "Finished parsing depletion library." << std::endl;
      if (_mico_xs_file_name != "")
      {
        _console << "Parsing microscopic cross-section library..." << std::endl;
        loadOpenMCMicoXSXML();
        _console << "Finished microscopic cross-section library." << std::endl;
      }
      else
        _console
            << COLOR_YELLOW
            << "Microscopic cross-sections have not been provided. Depletion calculations will "
               "only support radioactive decay.\n"
            << COLOR_DEFAULT;
      break;
    default:
      mooseError("Unsupported file source.");
      break;
  }

  _console << "-----------------------------------------------------" << std::endl;
}

void
DepletionLibraryAction::act()
{
  if (_current_task == "add_user_object" && getParam<bool>("add_data_userobject"))
  {
    auto params = _factory.getValidParams("DepletionDataProvider");
    params.set<std::string>("depletion_file_name") = getParam<std::string>("depletion_file");
    params.set<std::string>("xs_file_name") = getParam<std::string>("cross_section_file");
    params.set<MooseEnum>("depletion_file_source") = getParam<MooseEnum>("depletion_file_source");
    params.set<std::vector<Real>>("group_boundaries") = _group_bounds;

    _problem->addUserObject(
        "DepletionDataProvider", getParam<std::string>("depletion_uo_name"), params);
  }
}

void
DepletionLibraryAction::getNuclidesFromInitial(std::vector<std::string> & nuclides) const
{
  using namespace NuclearData;

  std::queue<std::string> traverse;
  for (const auto & nuclide : nuclides)
    traverse.push(nuclide);

  while (!traverse.empty())
  {
    if (_nuclide_list.count(traverse.front()) > 0u)
    {
      const auto nuclide = _nuclide_list.at(traverse.front());
      for (const auto & decay : nuclide.getDecays())
      {
        if (std::find(nuclides.begin(), nuclides.end(), decay._target) == nuclides.end())
        {
          if (decay._target != "")
          {
            traverse.push(decay._target);
            nuclides.emplace_back(decay._target);
          }
        }
      }

      for (const auto & rxn : nuclide.getReactions())
      {
        if (std::find(nuclides.begin(), nuclides.end(), rxn._target) == nuclides.end())
        {
          if (rxn._target != "")
          {
            traverse.push(rxn._target);
            nuclides.emplace_back(rxn._target);
          }
        }
      }
    }

    traverse.pop();
  }
}

void
DepletionLibraryAction::loadOpenMCDepletionXML()
{
  using namespace NuclearData;

  pugi::xml_document doc;
  if (!doc.load_file(_depletion_file_name.c_str()))
    mooseError("Failed to load the following depletion file: " + _depletion_file_name);

  // Traverse through the document and fetch the nuclides.
  std::vector<Real> storage;
  for (auto & nuclide_node : doc.child("depletion_chain"))
  {
    // Name and half-life.
    const std::string name(nuclide_node.attribute("name").as_string());
    Real half_life = nuclide_node.attribute("half_life").as_double(-1.0);

    _nuclide_list.emplace(name, Nuclide(name, half_life));

    for (auto & data_node : nuclide_node)
    {
      // Parse decay data.
      if (std::string(data_node.name()) == "decay")
      {
        if (std::string(data_node.attribute("type").as_string()) == "it" ||
            std::string(data_node.attribute("type").as_string()) == "IT")
        {
          const std::string target(data_node.attribute("target").as_string());
          _nuclide_list.at(name).addDecay(
              Decay::Mode::IT, data_node.attribute("branching_ratio").as_double(), target);

          if (_nuclide_decay_parents.count(target) == 0u)
            _nuclide_decay_parents.emplace(target, std::vector<std::string>());

          _nuclide_decay_parents.at(target).emplace_back(name);
        }
        else if (std::string(data_node.attribute("type").as_string()) == "beta-")
        {
          const std::string target(data_node.attribute("target").as_string());
          _nuclide_list.at(name).addDecay(
              Decay::Mode::BetaMinus, data_node.attribute("branching_ratio").as_double(), target);

          if (_nuclide_decay_parents.count(target) == 0u)
            _nuclide_decay_parents.emplace(target, std::vector<std::string>());

          _nuclide_decay_parents.at(target).emplace_back(name);
        }
        else if (std::string(data_node.attribute("type").as_string()) == "beta+")
        {
          const std::string target(data_node.attribute("target").as_string());
          _nuclide_list.at(name).addDecay(
              Decay::Mode::BetaPlus, data_node.attribute("branching_ratio").as_double(), target);

          if (_nuclide_decay_parents.count(target) == 0u)
            _nuclide_decay_parents.emplace(target, std::vector<std::string>());

          _nuclide_decay_parents.at(target).emplace_back(name);
        }
        else if (std::string(data_node.attribute("type").as_string()) == "ec/beta+")
        {
          const std::string target(data_node.attribute("target").as_string());
          _nuclide_list.at(name).addDecay(
              Decay::Mode::ECBetaPlus, data_node.attribute("branching_ratio").as_double(), target);

          if (_nuclide_decay_parents.count(target) == 0u)
            _nuclide_decay_parents.emplace(target, std::vector<std::string>());

          _nuclide_decay_parents.at(target).emplace_back(name);
        }
        else if (std::string(data_node.attribute("type").as_string()) == "alpha")
        {
          const std::string target(data_node.attribute("target").as_string());
          _nuclide_list.at(name).addDecay(
              Decay::Mode::Alpha, data_node.attribute("branching_ratio").as_double(), target);

          if (_nuclide_decay_parents.count(target) == 0u)
            _nuclide_decay_parents.emplace(target, std::vector<std::string>());

          _nuclide_decay_parents.at(target).emplace_back(name);
        }
        else if (std::string(data_node.attribute("type").as_string()) == "sf" ||
                 std::string(data_node.attribute("type").as_string()) == "SF")
        {
          const std::string target(data_node.attribute("target").as_string());
          _nuclide_list.at(name).addDecay(
              Decay::Mode::SF, data_node.attribute("branching_ratio").as_double(), target);

          if (_nuclide_decay_parents.count(target) == 0u)
            _nuclide_decay_parents.emplace(target, std::vector<std::string>());

          _nuclide_decay_parents.at(target).emplace_back(name);
        }
        else if (std::string(data_node.attribute("type").as_string()) == "beta-,alpha")
        {
          const std::string target(data_node.attribute("target").as_string());
          _nuclide_list.at(name).addDecay(Decay::Mode::BetaMinusAlpha,
                                          data_node.attribute("branching_ratio").as_double(),
                                          target);

          if (_nuclide_decay_parents.count(target) == 0u)
            _nuclide_decay_parents.emplace(target, std::vector<std::string>());

          _nuclide_decay_parents.at(target).emplace_back(name);
        }
        else if (std::string(data_node.attribute("type").as_string()) == "beta-,n")
        {
          const std::string target(data_node.attribute("target").as_string());
          _nuclide_list.at(name).addDecay(Decay::Mode::BetaMinusNeutron,
                                          data_node.attribute("branching_ratio").as_double(),
                                          target);

          if (_nuclide_decay_parents.count(target) == 0u)
            _nuclide_decay_parents.emplace(target, std::vector<std::string>());

          _nuclide_decay_parents.at(target).emplace_back(name);
        }
        else if (_warnings)
          _console << COLOR_YELLOW << "Unsupported decay for " << name << ": '"
                   << data_node.attribute("type").as_string() << "'.\n"
                   << COLOR_DEFAULT;
      }

      // Parse particle sources (per decay).
      if (std::string(data_node.name()) == "source")
      {
        if (std::string(data_node.attribute("type").as_string()) == "discrete")
        {
          if (std::string(data_node.attribute("particle").as_string()) == "photon")
          {
            parseToVector(std::string(data_node.child("parameters").child_value()), storage);
            if (storage.size() > 0u)
              _nuclide_list.at(name).addSource(Particletype::GammaPhoton, storage);
            storage.clear();
          }
          else if (std::string(data_node.attribute("particle").as_string()) == "neutron")
          {
            parseToVector(std::string(data_node.child("parameters").child_value()), storage);
            if (storage.size() > 0u)
              _nuclide_list.at(name).addSource(Particletype::Neutron, storage);
            storage.clear();
          }
          else if (_warnings)
            _console << "Source particle '" << data_node.attribute("particle").as_string()
                     << "' is not supported.\n";
        }
        else if (_warnings)
          _console << COLOR_YELLOW << "Source distribution '"
                   << data_node.attribute("type").as_string() << "' is not supported.\n"
                   << COLOR_DEFAULT;
      }

      // Parsing reaction data.
      if (std::string(data_node.name()) == "reaction")
      {
        if (std::string(data_node.attribute("type").as_string()) == "(n,gamma)")
        {
          const auto target = _nuclide_list.at(name).addReaction(
              Reaction::Mode::NGamma,
              data_node.attribute("branching_ratio").as_double(1.0),
              std::string(data_node.attribute("target").as_string()),
              0,
              1,
              data_node.attribute("Q").as_double());

          if (_nuclide_activation_parents.count(target) == 0u)
            _nuclide_activation_parents.emplace(target, std::vector<std::string>());

          _nuclide_activation_parents.at(target).emplace_back(name);
        }
        else if (std::string(data_node.attribute("type").as_string()) == "(n,p)")
        {
          const auto target = _nuclide_list.at(name).addReaction(
              Reaction::Mode::NProton,
              data_node.attribute("branching_ratio").as_double(1.0),
              std::string(data_node.attribute("target").as_string()),
              -1,
              0,
              data_node.attribute("Q").as_double());

          if (_nuclide_activation_parents.count(target) == 0u)
            _nuclide_activation_parents.emplace(target, std::vector<std::string>());

          _nuclide_activation_parents.at(target).emplace_back(name);
        }
        else if (std::string(data_node.attribute("type").as_string()) == "(n,a)")
        {
          const auto target = _nuclide_list.at(name).addReaction(
              Reaction::Mode::NAlpha,
              data_node.attribute("branching_ratio").as_double(1.0),
              std::string(data_node.attribute("target").as_string()),
              -2,
              -1,
              data_node.attribute("Q").as_double());

          if (_nuclide_activation_parents.count(target) == 0u)
            _nuclide_activation_parents.emplace(target, std::vector<std::string>());

          _nuclide_activation_parents.at(target).emplace_back(name);
        }
        else if (std::string(data_node.attribute("type").as_string()) == "(n,2n)")
        {
          const auto target = _nuclide_list.at(name).addReaction(
              Reaction::Mode::N2N,
              data_node.attribute("branching_ratio").as_double(1.0),
              std::string(data_node.attribute("target").as_string()),
              0,
              -1,
              data_node.attribute("Q").as_double());

          if (_nuclide_activation_parents.count(target) == 0u)
            _nuclide_activation_parents.emplace(target, std::vector<std::string>());

          _nuclide_activation_parents.at(target).emplace_back(name);
        }
        else if (std::string(data_node.attribute("type").as_string()) == "(n,3n)")
        {
          const auto target = _nuclide_list.at(name).addReaction(
              Reaction::Mode::N3N,
              data_node.attribute("branching_ratio").as_double(1.0),
              std::string(data_node.attribute("target").as_string()),
              0,
              -2,
              data_node.attribute("Q").as_double());

          if (_nuclide_activation_parents.count(target) == 0u)
            _nuclide_activation_parents.emplace(target, std::vector<std::string>());

          _nuclide_activation_parents.at(target).emplace_back(name);
        }
        else if (std::string(data_node.attribute("type").as_string()) == "(n,4n)")
        {
          const auto target = _nuclide_list.at(name).addReaction(
              Reaction::Mode::N4N,
              data_node.attribute("branching_ratio").as_double(1.0),
              std::string(data_node.attribute("target").as_string()),
              0,
              -3,
              data_node.attribute("Q").as_double());

          if (_nuclide_activation_parents.count(target) == 0u)
            _nuclide_activation_parents.emplace(target, std::vector<std::string>());

          _nuclide_activation_parents.at(target).emplace_back(name);
        }
        else if (_warnings)
          _console << COLOR_YELLOW << "Unsupported reaction for " << name << ": '"
                   << data_node.attribute("type").as_string() << "'.\n"
                   << COLOR_DEFAULT;
      }
    }
  }
}

void
DepletionLibraryAction::loadOpenMCMicoXSXML()
{
  using namespace NuclearData;

  pugi::xml_document doc;
  if (!doc.load_file(_mico_xs_file_name.c_str()))
    mooseError("Failed to load the following microscopic cross-section file: " +
               _mico_xs_file_name);

  auto chain = doc.child("depletion_chain");
  if (std::string(chain.attribute("energy_units").as_string()) == "eV")
    _convert_to_eV = 1.0;
  else if (std::string(chain.attribute("energy_units").as_string()) == "keV")
    _convert_to_eV = 1.0e3;
  else if (std::string(chain.attribute("energy_units").as_string()) == "MeV")
    _convert_to_eV = 1.0e6;
  else
    _console << COLOR_YELLOW
             << "Unsupported units of energy: " << chain.attribute("energy_units").as_string()
             << ". Defauling to eV.\n"
             << COLOR_DEFAULT;

  if (std::string(chain.attribute("xs_units").as_string()) == "barns")
    _xs_units = XSUnits::Barns;
  else if (std::string(chain.attribute("xs_units").as_string()) == "cm^-2")
    _xs_units = XSUnits::InvCm2;
  else
    _console << COLOR_YELLOW << "Cross-section units cannot be detected. Defauling to barns.\n"
             << COLOR_DEFAULT;

  // Parse the number of groups and group boundaries.
  _num_groups = chain.attribute("num_groups").as_uint(0u);
  parseToVector(std::string(chain.attribute("group_bounds").as_string()), _group_bounds);
  for (auto & bnd : _group_bounds)
    bnd *= _convert_to_eV;

  if (_num_groups + 1 != _group_bounds.size())
    mooseError("The number of groups provided is not consistent with the group boundaries.");

  // Traverse through the document and fetch the nuclides.
  std::vector<Real> storage;
  for (auto & nuclide_node : doc.child("depletion_chain"))
  {
    // Nuclide name.
    const std::string name(nuclide_node.attribute("name").as_string());

    if (_nuclide_list.count(name) == 0u)
    {
      if (_warnings)
        _console
            << COLOR_YELLOW << "Nuclide " << name
            << " does not exist in the depletion chain file. This nuclide will be skipped when "
               "parsing cross-sections.\n"
            << COLOR_DEFAULT;
      continue;
    }

    for (auto & data_node : nuclide_node)
    {
      // Parsing reaction data.
      if (std::string(data_node.name()) == "reaction")
      {
        if (std::string(data_node.attribute("type").as_string()) == "(n,gamma)")
        {
          parseToVector(std::string(data_node.attribute("mgxs").as_string()), storage);
          if (_xs_units == XSUnits::Barns)
            for (auto & xs : storage)
              xs *= 1e-24;

          if (!_nuclide_list.at(name).addReactionCrossSections(Reaction::Mode::NGamma, storage) &&
              _warnings)
            _console << COLOR_YELLOW << "Failed to add cross-sections to " << name
                     << " for the reaction " << std::string(data_node.attribute("type").as_string())
                     << COLOR_DEFAULT << "\n";
          storage.clear();
        }
        else if (std::string(data_node.attribute("type").as_string()) == "(n,p)")
        {
          parseToVector(std::string(data_node.attribute("mgxs").as_string()), storage);
          if (_xs_units == XSUnits::Barns)
            for (auto & xs : storage)
              xs *= 1e-24;

          if (!_nuclide_list.at(name).addReactionCrossSections(Reaction::Mode::NProton, storage) &&
              _warnings)
            _console << COLOR_YELLOW << "Failed to add cross-sections to " << name
                     << " for the reaction " << std::string(data_node.attribute("type").as_string())
                     << COLOR_DEFAULT << "\n";
          storage.clear();
        }
        else if (std::string(data_node.attribute("type").as_string()) == "(n,a)")
        {
          parseToVector(std::string(data_node.attribute("mgxs").as_string()), storage);
          if (_xs_units == XSUnits::Barns)
            for (auto & xs : storage)
              xs *= 1e-24;

          if (!_nuclide_list.at(name).addReactionCrossSections(Reaction::Mode::NAlpha, storage) &&
              _warnings)
            _console << COLOR_YELLOW << "Failed to add cross-sections to " << name
                     << " for the reaction " << std::string(data_node.attribute("type").as_string())
                     << COLOR_DEFAULT << "\n";
          storage.clear();
        }
        else if (std::string(data_node.attribute("type").as_string()) == "(n,2n)")
        {
          parseToVector(std::string(data_node.attribute("mgxs").as_string()), storage);
          if (_xs_units == XSUnits::Barns)
            for (auto & xs : storage)
              xs *= 1e-24;

          if (!_nuclide_list.at(name).addReactionCrossSections(Reaction::Mode::N2N, storage) &&
              _warnings)
            _console << COLOR_YELLOW << "Failed to add cross-sections to " << name
                     << " for the reaction " << std::string(data_node.attribute("type").as_string())
                     << COLOR_DEFAULT << "\n";
          storage.clear();
        }
        else if (std::string(data_node.attribute("type").as_string()) == "(n,3n)")
        {
          parseToVector(std::string(data_node.attribute("mgxs").as_string()), storage);
          if (_xs_units == XSUnits::Barns)
            for (auto & xs : storage)
              xs *= 1e-24;

          if (!_nuclide_list.at(name).addReactionCrossSections(Reaction::Mode::N3N, storage) &&
              _warnings)
            _console << COLOR_YELLOW << "Failed to add cross-sections to " << name
                     << " for the reaction " << std::string(data_node.attribute("type").as_string())
                     << COLOR_DEFAULT << "\n";
          storage.clear();
        }
        else if (std::string(data_node.attribute("type").as_string()) == "(n,4n)")
        {
          parseToVector(std::string(data_node.attribute("mgxs").as_string()), storage);
          if (_xs_units == XSUnits::Barns)
            for (auto & xs : storage)
              xs *= 1e-24;

          if (!_nuclide_list.at(name).addReactionCrossSections(Reaction::Mode::N4N, storage) &&
              _warnings)
            _console << COLOR_YELLOW << "Failed to add cross-sections to " << name
                     << " for the reaction " << std::string(data_node.attribute("type").as_string())
                     << COLOR_DEFAULT << "\n";
          storage.clear();
        }
        else
          _console << "Unsupported reaction for " << name << ": '"
                   << data_node.attribute("type").as_string() << "'.\n";
      }
    }
  }
}

void
DepletionLibraryAction::parseToVector(const std::string & string_rep, std::vector<Real> & real_rep)
{
  if (string_rep.size() == 0u)
    return;

  auto current_delim_pos = string_rep.find(" ");
  std::size_t offset = 0u;
  auto previous_delim_pos = 0u;
  while (current_delim_pos != std::string::npos)
  {
    real_rep.emplace_back(std::stod(string_rep.substr(
        previous_delim_pos + offset, current_delim_pos - (previous_delim_pos + offset))));
    previous_delim_pos = current_delim_pos;
    current_delim_pos = string_rep.find(" ", current_delim_pos + 1);
    offset = 1u;
  }

  real_rep.emplace_back(std::stod(string_rep.substr(previous_delim_pos)));
}
