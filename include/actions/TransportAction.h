#pragma once

#include "GnatBaseAction.h"

#include <unordered_map>
#include <vector>

// TODO:
// - This action + syntax is in desperate need of a refactor.
// - Finish function initial conditions.
class TransportAction : public GnatBaseAction
{
public:
  static InputParameters validParams();

  TransportAction(const InputParameters & params);

  virtual void act() override;

protected:
  static bool vecEquals(const RealVectorValue & first,
                        const RealVectorValue & second,
                        const Real & tol = libMesh::TOLERANCE);

  // Helper member function to initialize SN quadrature parameters.
  void applyQuadratureParameters(InputParameters & params);

  // Member function to initialize common scheme parameters.
  void actCommon();

  // Individual act functions for each scheme.
  void actSAAFCFEM();
  void actDiffusion();
  void actTransfer();

  // Act function for setting up the multi-app provided uncollided flux treatment.
  void actUncollided();

  // Member functions to initialize the MOOSE objects required for all schemes.
  void modifyOutputs();
  void addSNUserObjects();
  void addSNBCs(const std::string & var_name, unsigned int g, unsigned int n);
  void addSNICs(const std::string & var_name, unsigned int g);
  void addAuxVariables(const std::string & var_name);
  // is_output = true as a default parameter.
  // TODO: Fix this for scattering moments.
  void addAuxKernels(const std::string & var_name, unsigned int g, unsigned int l, int m);

  // Member functions to initialize the MOOSE objects required for the
  // CGFEM-SAAF scheme.
  void addSAAFKernels(const std::string & var_name, unsigned int g, unsigned int n);
  void addSAAFDiracKernels(const std::string & var_name, unsigned int g, unsigned int n);

  // Member functions to initialize the MOOSE objects required for the
  // diffusion approximation scheme.
  void addDiffusionBCs(const std::string & var_name);
  void addDiffusionICs(const std::string & var_name, unsigned int g);
  void addDiffusionKernels(const std::string & var_name, unsigned int g);
  void addDiffusionDiracKernels(const std::string & var_name, unsigned int g);

  // Member functions to add transfers.
  void addTransfers(const std::string & to_var_name, const std::string & source_var_name);
  void addUncTransfers(const std::string & to_var_name, const std::string & source_var_name);

  // Member function to add conservative transfer post-processors.
  void addDestinationConservativePP(const std::string & to_var_name);
  void addSourceConservativePP(const std::string & source_var_name);

  const TransportScheme _transport_scheme;
  const Particletype _particle;
  const bool _is_eigen;

  // Number of discrete ordinates, flux moments, and quadrature parameters.
  const unsigned int _n_l;
  const unsigned int _n_c;
  unsigned int _num_flux_ordinates;
  unsigned int _num_group_moments;

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
  const std::string & _flux_moment_name;

  // Base names of the stored directional angular fluxes.
  // Phi_{g, l, m} = _angular_flux_name_g_l_m.
  const std::string & _angular_flux_name;

  // Boundary side-sets.
  const std::vector<BoundaryName> _vacuum_side_sets;
  const std::vector<BoundaryName> _source_side_sets;
  const std::vector<BoundaryName> _current_side_sets;
  const std::vector<BoundaryName> _reflective_side_sets;

  // Point source moments.
  const std::vector<Point> & _point_source_locations;
  std::vector<std::vector<Real>> _point_source_moments;
  const std::vector<unsigned int> & _point_source_anisotropy;

  // Boundary source properties.
  std::vector<std::vector<Real>> _boundary_source_moments;
  const std::vector<unsigned int> & _boundary_source_anisotropy;

  // Boundary current properties.
  std::vector<std::vector<Real>> _boundary_currents;
  const std::vector<unsigned int> & _boundary_current_anisotropy;

  // Volumetric sources.
  const std::vector<SubdomainName> & _volumetric_source_blocks;
  std::vector<std::vector<Real>> _volumetric_source_moments;
  const std::vector<unsigned int> & _volumetric_source_anisotropy;

  // Field sources (radionuclides in plumes and similar sources).
  const std::vector<SubdomainName> & _field_source_blocks;
  const std::vector<std::vector<VariableName>> & _field_source_moments;
  const std::vector<unsigned int> & _field_source_anisotropy;
  const std::vector<Real> & _field_source_scaling;

  // Multi-app properties.
  const MultiAppName & _from_multi_app_name;
  const std::vector<SubdomainName> & _from_subdomain_ids;
  const std::string & _source_flux_moment_names;

  // Uncollided flux multi-app properties.
  const MultiAppName & _uncollided_from_multi_app_name;
  const std::vector<SubdomainName> & _uncollided_from_subdomain_ids;
  const std::string & _uncollided_source_flux_moment_names;
  const bool _using_uncollided;

  // List of the names for all angular flux variables and flux moments.
  std::unordered_map<unsigned int, std::vector<VariableName>> _group_angular_fluxes;
  std::unordered_map<unsigned int, std::vector<VariableName>> _group_flux_moments;

  // Source scaling.
  Real _source_scale_factor;

  bool _var_init;
}; // class TransportAction
