#pragma once

#include "GnatBaseAction.h"

#include <unordered_map>
#include <vector>

// TODO:
// - Finish function and file initial conditions.
// - Fix negative fluxes in the SAAF implementation.
//   - Tempoary fix through maxing the angular flux while computing flux moments.
class NeutronTransportAction : public GnatBaseAction
{
public:
  static InputParameters validParams();

  NeutronTransportAction(const InputParameters & params);

  virtual void act() override;

protected:
  using Action::addRelationshipManagers;
  void addRelationshipManagers(Moose::RelationshipManagerType when_type) override;

  // Helper member function to initialize SN quadrature parameters.
  void applyQuadratureParameters(InputParameters & params);

  // Member function to initialize common scheme parameters.
  void initializeCommon();

  // Individual act functions for each scheme.
  void actSAAFCFEM();
  void actUpwindDFEM();
  void actDiffusion();

  // Member functions to initialize the MOOSE objects required for all schemes.
  void addOutputs();
  void addSNBCs(const std::string & var_name, unsigned int g, unsigned int n);
  void addSNICs(const std::string & var_name, unsigned int g);
  void addAuxVariables(const std::string & var_name);
  // is_output = true as a default parameter.
  // TODO: Fix this for scattering moments.
  void addAuxKernels(const std::string & var_name, unsigned int g, unsigned int l, int m);

  // Member functions to initialize the MOOSE objects required for the
  // CGFEM-SAAF scheme.
  void addSAAFKernels(const std::string & var_name, unsigned int g, unsigned int n);
  void addSAAFDiracKernels();

  // Member functions to initialize the MOOSE objects required for the
  // DGFEM-upwinding scheme.
  void addDGFEMKernels(const std::string & var_name, unsigned int g, unsigned int n);
  void addDGFEMDGKernels(const std::string & var_name, unsigned int n);
  void addDGFEMDiracKernels();

  // Member functions to initialize the MOOSE objects required for the
  // diffusion approximation scheme.
  void addDiffusionBCs(const std::string & var_name);
  void addDiffusionICs(const std::string & var_name, unsigned int g);
  void addDiffusionKernels(const std::string & var_name, unsigned int g);
  void addDiffusionDiracKernels();

  const Scheme _transport_scheme;

  // Number of discrete ordinates.
  const unsigned int _n_l;
  const unsigned int _n_c;
  unsigned int _num_flux_ordinates;

  // Base names of the stored directional angular fluxes.
  // Phi_{g, l, m} = _angular_flux_name_g_l_m.
  const std::string _angular_flux_name;

  // Boundary side-sets.
  const std::vector<BoundaryName> _vacuum_side_sets;
  const std::vector<BoundaryName> _source_side_sets;
  const std::vector<BoundaryName> _reflective_side_sets;

  // List of the names for all angular flux variables and flux moments.
  std::unordered_map<unsigned int, std::vector<VariableName>> _group_angular_fluxes;
  bool _var_init;
}; // class NeutronTransportAction
