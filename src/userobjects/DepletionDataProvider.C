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
  params.addParam<MooseEnum>("depletion_file_source",
                             MooseEnum("openmc", "openmc"),
                             "The depletion file generator. More depletion file sources will be "
                             "supported at a later date.");

  return params;
}

DepletionDataProvider::DepletionDataProvider(const InputParameters & parameters)
  : ThreadedGeneralUserObject(parameters),
    _file_source(getParam<MooseEnum>("depletion_file_source").getEnum<DepletionFileSource>()),
    _depletion_file_name((std::filesystem::path(_app.getInputFileName()).parent_path() /
                          std::filesystem::path(getParam<std::string>("depletion_file_name")))
                             .string())
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
  if (!doc.load_file(_depletion_file_name.c_str()))
    mooseError("Failed to load the depletion file: " + _depletion_file_name);

  // Traverse through the document and fetch the nuclides.
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
        else
          _console << "Unknown decay '" << data_node.attribute("type").as_string() << "'.\n";
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
