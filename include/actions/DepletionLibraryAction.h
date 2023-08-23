#pragma once

#include "Action.h"

#include "NuclearData.h"
#include "Nuclide.h"

class DepletionLibraryAction : public Action
{
public:
  static InputParameters validParams();

  DepletionLibraryAction(const InputParameters & parameters);

  virtual void act() final {}

  bool hasNuclide(const std::string & nuclide) const { return _nuclide_list.count(nuclide) > 0; }

  const NuclearData::Nuclide & getNuclide(const std::string & nuclide) const
  {
    if (_nuclide_list.count(nuclide) == 0)
      mooseError("Depletion data does not exist for '" + nuclide + "'.");

    return _nuclide_list.at(nuclide);
  }
  const std::vector<std::string> & getDecayParents(const std::string & nuclide) const
  {
    if (_nuclide_decay_parents.count(nuclide) == 0u)
      return _no_children;
    return _nuclide_decay_parents.at(nuclide);
  }
  const std::vector<std::string> & getActivationParents(const std::string & nuclide) const
  {
    if (_nuclide_activation_parents.count(nuclide) == 0u)
      return _no_children;
    return _nuclide_activation_parents.at(nuclide);
  }

  // Get all nuclides in the depletion system given a list of initial nuclides (in-place).
  void getNuclidesFromInitial(std::vector<std::string> & nuclides) const;
  // Get the number of energy groups that this depletion library supports.
  unsigned int numGroups() const { return _num_groups; }

protected:
  void loadOpenMCDepletionXML();
  void loadOpenMCMicoXSXML();

  void parseToVector(const std::string & string_rep, std::vector<Real> & real_rep);

  const enum class DepletionFileSource { OpenMC = 0u } _file_source;

  const std::string _depletion_file_name;
  const std::string _mico_xs_file_name;

  const bool _warnings;

  std::unordered_map<std::string, NuclearData::Nuclide> _nuclide_list;

  // Pre-process the various decay and activation parents for each nuclide.
  std::unordered_map<std::string, std::vector<std::string>> _nuclide_decay_parents;
  std::unordered_map<std::string, std::vector<std::string>> _nuclide_activation_parents;
  const std::vector<std::string> _no_children;

  // Cross-section information.
  enum class XSUnits
  {
    InvCm2 = 0u,
    Barns = 1u
  } _xs_units;
  enum class EnergyUnits
  {
    eV = 0u,
    keV = 1u,
    MeV = 2u
  } _energy_units;

  unsigned int _num_groups;
  std::vector<Real> _group_bounds;
}; // class DepletionLibraryAction
