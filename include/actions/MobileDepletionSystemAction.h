#pragma once

#include "GnatBase.h"
#include "GnatBaseAction.h"

class NSFVAction;
class DepletionLibraryAction;

// A class which sets up depletion problems using a mass/number fraction based approach. More
// effective for problems where a fluid undergoes neutron transmutation and decay.
class MobileDepletionSystemAction : public GnatBaseAction
{
public:
  static InputParameters validParams();

  MobileDepletionSystemAction(const InputParameters & parameters);

  virtual void act() override;

protected:
  // Add elements using their natural isotopic abundance.
  void addElementsAndNuclides();
  // Build the depletion system.
  void buildDepletionSystem();
  // Fetch properties from the coupled TransportSystem.
  void fetchTransportProperties();

  // Modify outputs to remove the mass fraction primal variables.
  void modifyOutputs();

  void applyIsotopeParameters(InputParameters & params, bool apply_density = true);
  void addICs(const std::string & nuclide_var_name);
  void addMaterials(const std::string & nuclide_var_name);

  // Finite element variables.
  void addKernels(const std::string & nuclide_var_name);
  void addBCs(const std::string & nuclide_var_name);

  // Finite volume variables.
  void addFVVariables(const std::string & nuclide_var_name);
  void addFVKernels(const std::string & nuclide_var_name);
  void addFVBCs(const std::string & nuclide_var_name);

  // Auxvariables.
  void addAuxVariables(const std::string & nuclide_var_name);
  void addAuxKernels(const std::string & nuclide_var_name);

  // Auxvariables and auxkernels for radiation transport source terms.
  void addRadiationAuxVariables();
  void addRadiationAuxKernels();

  unsigned int _mesh_dims;
  bool _using_moose_ns;

  // The nuclide scheme.
  enum class NuclideScheme
  {
    SUPGFE = 0u,
    FV = 1u
  } _scheme;

  // The coupled Navier-Stokes finite volume action.
  const NSFVAction * _coupled_ns_fv;

  // The coupled Depletion Library.
  const DepletionLibraryAction * _coupled_depletion_lib;

  // The coupled TransportSystem properties.
  const std::string & _transport_system;
  bool _has_transport_system;
  unsigned int _num_groups;
  std::vector<VariableName> _group_flux_moments;

  bool _first_action;

  // The nuclides in the system. Name is paired with the associated density.
  const std::vector<std::string> & _elements;
  const std::vector<Real> & _element_atom_fractions;
  const std::vector<std::string> & _extra_nuclides;
  const std::vector<Real> & _extra_nuclide_atom_fractions;
  std::unordered_map<std::string, Real> _total_nuclide_list;

  // The boundary conditions in the system.
  const std::vector<BoundaryName> & _inlet_boundaries;
  const std::vector<std::vector<Real>> & _inlet_atom_fractions;
  const std::vector<BoundaryName> & _outlet_boundaries;
  std::vector<std::unordered_map<std::string, Real>> _boundary_mass_fractions;

  // For particle decay sources.
  const std::string & _photon_source_prefix;
  const std::string & _neutron_source_prefix;
  const std::vector<Real> & _photon_group_boundaries;
  const std::vector<Real> & _neutron_group_boundaries;

  // Debug options.
  const std::vector<std::string> & _debug_filter_nuclides;
}; // class MobileDepletionSystemAction
