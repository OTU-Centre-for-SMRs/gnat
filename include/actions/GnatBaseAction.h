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
  void debugOutput(const std::string & level0 = "",
                   const std::string & level1 = "");

  // The execution type (steady-state or transient) and the debug output verbosity.
  ExecutionType _exec_type;
  DebugVerbosity _debug_level;

  // Number of spectral energy groups.
  unsigned int _num_groups;
  // Maximum degree of the flux moments requested. By default the transport
  // action prepares the scalar flux and that's it.
  // The Gauss-Legendre-Chebyshev quadrature with half-circle orders
  // n_{g} and n_{c} can fully integrate a spherical harmonics expansion of
  // degree N IF: n_{g} <= N + 1 AND n_{c} <= N + 1.
  unsigned int _max_eval_anisotropy;

  // Base names of the stored angular flux moments.
  // Psi_{g, n} = _flux_moment_name_g_n.
  std::string _flux_moment_name;

  // The coordinate system type and dimensionality.
  ProblemType _p_type;
  unsigned int _num_group_moments;

  // Blocks to apply the neutron activation simulations.
  std::set<SubdomainID> _subdomain_ids;

  // List of the names for all flux moments.
  std::unordered_map<unsigned int, std::vector<VariableName>> _group_flux_moments;
};
