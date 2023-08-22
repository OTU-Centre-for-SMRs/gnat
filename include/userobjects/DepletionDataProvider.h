#pragma once

#include "ThreadedGeneralUserObject.h"

#include "NuclearData.h"
#include "Nuclide.h"

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

protected:
  void loadOpenMCXML();

  const enum class DepletionFileSource { OpenMC = 0u } _file_source;

  const std::string _depletion_file_name;

  std::unordered_map<std::string, NuclearData::Nuclide> _nuclide_list;
}; // class DepletionDataProvider.
