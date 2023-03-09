#pragma once

#include "Action.h"

#include "GnatBase.h"

class GnatBaseAction : public Action
{
public:
  static InputParameters validParams();

  GnatBaseAction(const InputParameters & params);

protected:
  // Helper member function to initialize base parameters.
  void initializeBase();
  // Helper member function for debug output.
  void debugOutput(const std::string & level0 = "", const std::string & level1 = "");

  // Add a variables.
  void addVariable(const std::string & var_name);

  // The execution type (steady-state or transient) and the debug output verbosity.
  ExecutionType _exec_type;
  DebugVerbosity _debug_level;

  // The coordinate system type and dimensionality.
  ProblemType _p_type;

  // Blocks to apply the neutron activation simulations.
  std::set<SubdomainID> _subdomain_ids;
};
