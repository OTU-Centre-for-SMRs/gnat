#pragma once

#include "GnatBase.h"
#include "GnatBaseAction.h"

class AddFileIsotopesAction : public GnatBaseAction
{
public:
  static InputParameters validParams();

  AddFileIsotopesAction(const InputParameters & parameters);

  virtual void act() override;

protected:
  struct NuclideProperties
  {
    Real _diffusion;
    Real _half_life;
    std::vector<Real> _sigma_a_g;

    std::vector<Real> _parent_decay_constants;
    std::vector<Real> _parent_branching_fractions;
    std::vector<Real> _sigma_a_g_parent;

    std::vector<VariableName> _decay_parents;
    std::vector<VariableName> _activation_parents;

    NuclideProperties() : _diffusion(0.0), _half_life(0.0) {}
  };

  void applyIsotopeParameters(InputParameters & params);
  void addICs(const std::string & var_name);
  void addKernels(const std::string & var_name, const NuclideProperties & properties);
  void addMaterials(const std::string & var_name, const NuclideProperties & properties);

  void parseProperty(const PropertyType & type, const std::string & property_file);
  void parseGnatProperty(const PropertyType & type, const std::string & property_file);
  void parseOpenMCProperty(const PropertyType & type, const std::string & property_file);

  // Properties for individual isotopes.
  std::unordered_map<std::string, NuclideProperties> _nuclide_properties;

  CrossSectionSource _xs_source;
  const std::string & _xs_file_name;
  const std::string & _xs_source_material_id;

  HalfLifeUnits _hl_units;
  const std::string & _nuclide_prop_file_name;

  bool _first_action;
}; // class AddFileIsotopesAction
