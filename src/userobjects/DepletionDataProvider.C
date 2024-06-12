#include "DepletionDataProvider.h"

#include <filesystem>
#include "pugixml.h"

registerMooseObject("GnatApp", DepletionDataProvider);

InputParameters
DepletionDataProvider::validParams()
{
  auto params = ThreadedGeneralUserObject::validParams();
  params.addClassDescription("A user object which provides and maintains a depletion library for "
                             "both stationary depletion and mobile depletion.");
  params.addRequiredParam<std::string>("depletion_file_name", "The name of the depletion file.");
  params.addRequiredParam<std::string>("xs_file_name", "The name of the depletion xs file.");
  params.addParam<MooseEnum>("depletion_file_source",
                             MooseEnum("openmc", "openmc"),
                             "The depletion file generator. More depletion file sources will be "
                             "supported at a later date.");
  params.addRequiredParam<std::vector<Real>>("group_boundaries",
                                             "The group structure (including 0.0 eV)");
  params.addParam<bool>(
      "show_warnings", false, "Whether or not this object should show warnings or not.");

  return params;
}

DepletionDataProvider::DepletionDataProvider(const InputParameters & parameters)
  : ThreadedGeneralUserObject(parameters),
    _file_source(getParam<MooseEnum>("depletion_file_source").getEnum<DepletionFileSource>()),
    _depletion_file((std::filesystem::path(_app.getLastInputFileName()).parent_path() /
                     std::filesystem::path(getParam<std::string>("depletion_file_name")))
                        .string()),
    _xs_file((std::filesystem::path(_app.getLastInputFileName()).parent_path() /
              std::filesystem::path(getParam<std::string>("xs_file_name")))
                 .string()),
    _group_bounds(getParam<std::vector<Real>>("group_boundaries")),
    _warnings(getParam<bool>("show_warnings"))
{
  switch (_file_source)
  {
    case DepletionFileSource::OpenMC:
      loadOpenMCXML();
      break;
    default:
      mooseError("Unsupported file source.");
      break;
  }
}

void
DepletionDataProvider::loadOpenMCXML()
{
  using namespace NuclearData;

  pugi::xml_document doc;
  if (!doc.load_file(_depletion_file.c_str()))
    mooseError("Failed to load the depletion file: " + _depletion_file);

  // Traverse through the document and fetch the nuclides.
  std::vector<Real> storage;
  for (auto & nuclide_node : doc.child("depletion_chain"))
  {
    // Name and half-life.
    std::string name(nuclide_node.attribute("name").as_string());
    Real half_life = nuclide_node.attribute("half_life").as_double(-1.0);

    _nuclide_list.emplace(name, Nuclide(name, half_life));

    for (auto & data_node : nuclide_node)
    {
      // Parse decay data.
      if (std::string(data_node.name()) == "decay")
      {
        if (std::string(data_node.attribute("type").as_string()) == "it" ||
            std::string(data_node.attribute("type").as_string()) == "IT")
          _nuclide_list.at(name).addDecay(Decay::Mode::IT,
                                          data_node.attribute("branching_ratio").as_double(),
                                          std::string(data_node.attribute("target").as_string()));
        else if (std::string(data_node.attribute("type").as_string()) == "beta-")
          _nuclide_list.at(name).addDecay(Decay::Mode::BetaMinus,
                                          data_node.attribute("branching_ratio").as_double(),
                                          std::string(data_node.attribute("target").as_string()));
        else if (std::string(data_node.attribute("type").as_string()) == "beta+")
          _nuclide_list.at(name).addDecay(Decay::Mode::BetaPlus,
                                          data_node.attribute("branching_ratio").as_double(),
                                          std::string(data_node.attribute("target").as_string()));
        else if (std::string(data_node.attribute("type").as_string()) == "ec/beta+")
          _nuclide_list.at(name).addDecay(Decay::Mode::ECBetaPlus,
                                          data_node.attribute("branching_ratio").as_double(),
                                          std::string(data_node.attribute("target").as_string()));
        else if (std::string(data_node.attribute("type").as_string()) == "alpha")
          _nuclide_list.at(name).addDecay(Decay::Mode::Alpha,
                                          data_node.attribute("branching_ratio").as_double(),
                                          std::string(data_node.attribute("target").as_string()));
        else if (std::string(data_node.attribute("type").as_string()) == "sf" ||
                 std::string(data_node.attribute("type").as_string()) == "SF")
          _nuclide_list.at(name).addDecay(Decay::Mode::SF,
                                          data_node.attribute("branching_ratio").as_double(),
                                          std::string(data_node.attribute("target").as_string()));
        else if (_warnings)
          _console << COLOR_YELLOW << "Unknown decay '" << data_node.attribute("type").as_string()
                   << "'.\n"
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

      // Parsing reaction data. Radiative capture (n, \gamma) are the only supported reactions at
      // the moment.
      if (std::string(data_node.name()) == "reaction" &&
          std::string(data_node.attribute("type").as_string()) == "(n,gamma)")
      {
        _nuclide_list.at(name).addReaction(Reaction::Mode::NGamma,
                                           data_node.attribute("branching_ratio").as_double(1.0),
                                           std::string(data_node.attribute("target").as_string()),
                                           0,
                                           1,
                                           data_node.attribute("Q").as_double());
      }
    }
  }
}

void
DepletionDataProvider::parseToVector(const std::string & string_rep, std::vector<Real> & real_rep)
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
