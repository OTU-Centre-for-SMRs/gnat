#pragma once

#include "ThreadedGeneralUserObject.h"

#include "NuclearData.h"
#include "Nuclide.h"

// TODO: Radionuclide sources. Branching factor for a photon emission = prob / sum(prob).

// A class which constructs and maintains depletion chains using user-provided depletion libraries.
// Currently supports depletion xml files generated using OpenMC. Eventually aims to support the
// isoXML file format used by Griffin.
class DepletionDataProvider : public ThreadedGeneralUserObject
{
public:
  static InputParameters validParams();

  DepletionDataProvider(const InputParameters & parameters);

  virtual void execute() final {}
  virtual void initialize() final {}
  virtual void finalize() final {}

  virtual void threadJoin(const UserObject &) final {}
  virtual void subdomainSetup() final {}

  // Group boundaries.
  bool validGroupBounds(unsigned int num_groups) const
  {
    return _group_bounds.size() - 1u == num_groups;
  }
  void getGroupBounds(unsigned int g, Real & top, Real & bottom) const
  {
    top = _group_bounds[g];
    bottom = _group_bounds[g + 1];
  }

  // Nuclide depletion data.
  bool hasNuclide(const std::string & nuclide) const
  {
    if (nuclide.find('_') != std::string::npos)
      return _nuclide_list.count(nuclide.substr(0u, nuclide.find('_'))) > 0u;
    else
      return _nuclide_list.count(nuclide) > 0u;
  }
  const NuclearData::Nuclide & getNuclide(const std::string & nuclide) const
  {
    return _nuclide_list.at(nuclide);
  }

protected:
  void loadOpenMCXML();
  void parseToVector(const std::string & string_rep, std::vector<Real> & real_rep);

  const enum class DepletionFileSource { OpenMC = 0u } _file_source;

  const std::string _depletion_file;
  const std::string _xs_file;

  std::unordered_map<std::string, NuclearData::Nuclide> _nuclide_list;

  const std::vector<Real> _group_bounds;

  const bool _warnings;
}; // class DepletionDataProvider.
