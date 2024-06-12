#pragma once

#include "AuxKernel.h"

class DepletionDataProvider;

// A class to compute the decay source term (gamma or neutron) from radionuclides.
class ParticleDecaySource : public AuxKernel
{
public:
  static InputParameters validParams();

  ParticleDecaySource(const InputParameters & parameters);

protected:
  virtual Real computeValue() override;

  const Particletype _particle;

  const std::vector<Real> & _group_bounds;
  const unsigned int _num_groups;
  const unsigned int _group_index;

  // The concentration variables in either variable (FE) or functor (FV) form.
  const bool _is_fe;
  std::vector<const ADVariableValue *> _var_nuclide_concentrations;
  std::vector<const Moose::Functor<ADReal> *> _fun_nuclide_concentrations;

  std::vector<std::string> _nuclide_names;

  const DepletionDataProvider * _data_provider;

  const bool _is_mass_density;
}; // class ParticleDecaySource
