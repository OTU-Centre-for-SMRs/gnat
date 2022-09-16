#pragma once

#include "GnatBaseAction.h"

class AddMobileIsotopeAction : public GnatBaseAction
{
public:
  static InputParameters validParams();

  AddMobileIsotopeAction(const InputParameters & params);

  virtual void act() override;

protected:
  using Action::addRelationshipManagers;
  void addRelationshipManagers(Moose::RelationshipManagerType when_type) override;

  void applyIsotopeParameters(InputParameters & params);
  void addKernels();
  void addMaterials();

  enum class HalfLifeUnits
  {
    Seconds = 0u,
    Minutes = 1u,
    Hours = 2u,
    Days = 3u,
    Years = 4u
  } _hl_units;

  // Variable names.
  const VariableName _isotope_name;
  const std::vector<VariableName> _decay_parents;
  const std::vector<VariableName> _activation_parents;

  const std::vector<VariableName> & _master_isotope_list;

  // Base diffusion coefficient. TODO: Different diffusion coefficient correlations?
  const Real _diffusion_coefficient_base;
  // Vector of microscopic absorption cross-sections.
  const std::vector<Real> _sigma_a;
  // Half-life.
  Real _half_life;

  // Decay parent isotope properties.
  std::vector<Real> _parent_decay_constants;
  const std::vector<Real> _parent_branching_fractions;

  // Activation parent isotope properties.
  const std::vector<Real> _parent_sigma_act;

  bool _first_action;
}; // class AddMobileIsotopeAction
