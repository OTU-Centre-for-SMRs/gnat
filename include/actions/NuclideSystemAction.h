#pragma once

#include "GnatBase.h"
#include "GnatBaseAction.h"

class NuclideSystemAction : public GnatBaseAction
{
public:
  static InputParameters validParams();

  NuclideSystemAction(const InputParameters & parameters);

  virtual void act() override;

protected:
  struct NuclideProperties
  {
    Real _diffusion;
    Real _half_life;
    Real _weight_fraction;
    std::vector<Real> _sigma_a_g;

    std::vector<Real> _parent_decay_constants;
    std::vector<Real> _parent_branching_fractions;
    std::vector<Real> _sigma_a_g_parent;

    std::vector<VariableName> _decay_parents;
    std::vector<VariableName> _activation_parents;

    NuclideProperties() : _diffusion(0.0), _half_life(-1.0), _weight_fraction(0.0) {}
  };

  void applyIsotopeParameters(InputParameters & params);
  void addICs(const std::string & var_name, const NuclideProperties & properties);
  void addKernels(const std::string & var_name, const NuclideProperties & properties);
  void addMaterials(const std::string & var_name, const NuclideProperties & properties);

  // Functions to assist with cross-section parsing.
  void parseProperty(const PropertyType & type, const std::string & property_file);
  void parseGnatProperty(const PropertyType & type, const std::string & property_file);
  void parseOpenMCProperty(const PropertyType & type, const std::string & property_file);

  // Function to assist with parsing the nuclear data.
  void parseNuclideSystem();
  void parseCrossSections();

  // Fetch properties from the coupled TransportSystem.
  void fetchTransportProperties();

  // The coupled TransportSystem properties.
  const std::string & _transport_system;
  bool _has_transport_system;
  unsigned int _num_groups;
  std::vector<VariableName> _group_flux_moments;

  Real _density;

  // Properties for individual isotopes.
  std::unordered_map<std::string, NuclideProperties> _nuclide_properties;

  CrossSectionSource _xs_source;
  CrossSectionType _xs_type;
  const std::string & _xs_file_name;
  const std::string & _xs_source_material_id;

  HalfLifeUnits _hl_units;
  const std::string & _nuclide_prop_file_name;

  bool _first_action;
}; // class NuclideSystemAction
