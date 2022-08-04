#pragma once

#include "Action.h"

#include "GnatBase.h"

#include <unordered_map>
#include <vector>

// TODO:
// - Complete the reflective boundary condition implementation.
// - Finish function and file initial conditions.
// - Fix negative fluxes in the SAAF implementation.
//   - Tempoary fix through maxing the angular flux while computing flux moments.
// - Fix the DGFEM implementation. Issue probably arises in the upwinding DG
//   kernel.
class NeutronTransportAction : public Action
{
public:
  static InputParameters validParams();

  NeutronTransportAction(const InputParameters & params);

  virtual void act() override;

protected:
  using Action::addRelationshipManagers;
  void addRelationshipManagers(Moose::RelationshipManagerType when_type) override;

  // Helper member function for debug output.
  void debugOutput(const std::string & level0 = "",
                   const std::string & level1 = "");
  // Helper member function to initialize SN quadrature parameters.
  void applyQuadratureParameters(InputParameters & params);

  // Member function to initialize common scheme parameters.
  void initializeCommon();

  // Individual act functions for each scheme.
  void actSAAFCFEM();
  void actUpwindDFEM();

  // Member functions to initialize the MOOSE objects required for all schemes.
  void addOutputs();
  void addVariable(const std::string & var_name);
  void addBCs(const std::string & var_name, unsigned int g, unsigned int n);
  void addICs(const std::string & var_name, unsigned int g, unsigned int n);
  void addAuxVariables(const std::string & var_name);
  // is_output = true as a default parameter.
  // TODO: Fix this for scattering moments.
  void addAuxKernels(const std::string & var_name, unsigned int g,
                     unsigned int l, int m);

  // Member functions to initialize the MOOSE objects required for the
  // CGFEM-SAAF scheme.
  void addSAAFKernels(const std::string & var_name, unsigned int g, unsigned int n);
  void addSAAFDiracKernels();

  // Member functions to initialize the MOOSE objects required for the
  // DGFEM-upwinding scheme.
  void addDGFEMKernels(const std::string & var_name, unsigned int g, unsigned int n);
  void addDGFEMDGKernels(const std::string & var_name, unsigned int g, unsigned int n);
  void addDGFEMDiracKernels();

  const Scheme _transport_scheme;
  const ExecutionType _exec_type;
  const DebugVerbosity _debug_level;

  std::set<SubdomainID> _subdomain_ids;

  // Number of spectral energy groups and number of discrete ordinates.
  const unsigned int _num_groups;
  const unsigned int _n_l;
  const unsigned int _n_c;

  // Maximum degree of the flux moments requested. By default the transport
  // action prepares the scalar flux and that's it.
  // The Gauss-Legendre-Chebyshev quadrature with half-circle orders
  // n_{g} and n_{c} can fully integrate a spherical harmonics expansion of
  // degree N IF: n_{g} <= N + 1 AND n_{c} <= N + 1.
  const unsigned int _max_eval_anisotropy;

  // Base names of the stored directional angular fluxes and angular flux moments.
  // As an example for both:
  // Psi_{g, n} = _flux_moment_name_g_n.
  // Phi_{g, l, m} = _angular_flux_name_g_l_m.
  const std::string _angular_flux_name;
  const std::string _flux_moment_name;

  // The coordinate system type and dimensionality.
  ProblemType _p_type;
  unsigned int _num_group_moments;
  unsigned int _num_flux_ordinates;

  // Boundary side-sets.
  const std::vector<BoundaryName> _vacuum_side_sets;
  const std::vector<BoundaryName> _source_side_sets;
  const std::vector<BoundaryName> _reflective_side_sets;

  // List of the names for all angular flux variables and flux moments.
  std::unordered_map<unsigned int, std::vector<VariableName>> _group_angular_fluxes;
  std::unordered_map<unsigned int, std::vector<VariableName>> _group_flux_moments;
  bool _var_init;
}; // class NeutronTransportAction
