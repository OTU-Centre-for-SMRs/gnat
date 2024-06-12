#include "ParticleDecaySource.h"

#include "DepletionDataProvider.h"

registerMooseObject("GnatApp", ParticleDecaySource);

InputParameters
ParticleDecaySource::validParams()
{
  auto params = AuxKernel::validParams();
  params.addClassDescription("Auxkernel which computes a particle source from the radioactive "
                             "decay of a nuclide concentration field. Useful for coupled radiation "
                             "/ mass transport for radiological consequence assessment.");

  params.addRequiredParam<MooseEnum>(
      "particle_type", MooseEnum("neutron photon"), "The particle this source represents.");
  params.addRequiredParam<std::vector<Real>>("group_boundaries",
                                             "The group structure (including 0.0 eV)");
  params.addRequiredParam<unsigned int>("group_index", "The current particle energy group.");

  params.addRequiredParam<bool>("is_fe", "Whether the mass fraction is finite element or not.");
  params.addCoupledVar("nuclide_vars", "The radionuclide concentration.");
  params.addParam<std::vector<MooseFunctorName>>("nuclide_funs", "The radionuclide concentration.");

  params.addParam<UserObjectName>("data_lib_name",
                                  "DepletionDataProviderUO",
                                  "The name of the depletion data provider userobject.");

  params.addParam<bool>(
      "is_number_density",
      true,
      "Whether the provided nuclide densities are number densities or mass densities.");

  return params;
}

ParticleDecaySource::ParticleDecaySource(const InputParameters & parameters)
  : AuxKernel(parameters),
    _particle(getParam<MooseEnum>("particle_type").getEnum<Particletype>()),
    _group_bounds(getParam<std::vector<Real>>("group_boundaries")),
    _num_groups(_group_bounds.size() - 1u),
    _group_index(getParam<unsigned int>("group_index")),
    _is_fe(getParam<bool>("is_fe")),
    _data_provider(&getUserObject<DepletionDataProvider>("data_lib_name", true)),
    _is_mass_density(!getParam<bool>("is_number_density"))
{
  if (_is_fe)
  {
    auto num_coupled = coupledComponents("nuclide_vars");

    // Validate concentration names to make sure we have a match for depletion data.
    for (unsigned int i = 0u; i < num_coupled; ++i)
    {
      if (!_data_provider->hasNuclide(coupledName("nuclide_vars", i)))
        mooseError("The nuclide " + coupledName("nuclide_vars", i) +
                   " does not exist in the depletion system.");
    }

    for (unsigned int i = 0u; i < num_coupled; ++i)
    {
      _var_nuclide_concentrations.emplace_back(&adCoupledValue("nuclide_vars", i));
      _nuclide_names.emplace_back(coupledName("nuclide_vars", i));
    }
  }
  else
  {
    // Validate concentration names to make sure we have a match for depletion data.
    for (const auto & fun_name : getParam<std::vector<MooseFunctorName>>("nuclide_funs"))
    {
      if (!_data_provider->hasNuclide(fun_name))
        mooseError("The nuclide " + fun_name + " does not exist in the depletion system.");
    }

    for (const auto & fun_name : getParam<std::vector<MooseFunctorName>>("nuclide_funs"))
    {
      _fun_nuclide_concentrations.emplace_back(&getFunctor<ADReal>(fun_name));
      _nuclide_names.emplace_back(fun_name);
    }
  }
}

Real
ParticleDecaySource::computeValue()
{
  constexpr Real n_avogadro = 6.0221408e23;

  const Real bot_bnd = _group_bounds[_group_index + 1];
  const Real top_bnd = _group_bounds[_group_index];

  Real val = 0.0;
  Real nuclide_val = 0.0;
  if (_is_fe)
  {
    for (unsigned int i = 0u; i < _nuclide_names.size(); ++i)
    {
      // Sum up the photon-specific decay constants for this group.
      for (const auto & source : _data_provider->getNuclide(_nuclide_names[i]).getSources())
      {
        if (source._particle == _particle)
        {
          for (unsigned int j = 0u; j < source._p_energies.size(); ++j)
          {
            if (source._p_energies[j] >= bot_bnd && source._p_energies[j] < top_bnd)
              nuclide_val += source._p_decay_constants[j];
          }
        }
      }
      // Multiply by number density.
      if (_is_mass_density)
        val += nuclide_val * MetaPhysicL::raw_value((*_var_nuclide_concentrations[i])[_qp]) *
               (n_avogadro / NuclearData::Nuclide::getAtomicMass(_nuclide_names[i]));
      else
        val += nuclide_val * MetaPhysicL::raw_value((*_var_nuclide_concentrations[i])[_qp]);

      nuclide_val = 0.0;
    }
  }
  else
  {
    for (unsigned int i = 0u; i < _nuclide_names.size(); ++i)
    {
      // Sum up the photon-specific decay constants for this group.
      for (const auto & source : _data_provider->getNuclide(_nuclide_names[i]).getSources())
      {
        if (source._particle == _particle)
        {
          for (unsigned int j = 0u; j < source._p_energies.size(); ++j)
          {
            if (source._p_energies[j] >= bot_bnd && source._p_energies[j] < top_bnd)
              nuclide_val += source._p_decay_constants[j];
          }
        }
      }
      // Multiply by number density.
      if (_is_mass_density)
        val += nuclide_val *
               MetaPhysicL::raw_value(
                   (*_fun_nuclide_concentrations[i])(makeElemArg(_current_elem), 0u)) *
               (n_avogadro / NuclearData::Nuclide::getAtomicMass(_nuclide_names[i]));
      else
        val += nuclide_val * MetaPhysicL::raw_value(
                                 (*_fun_nuclide_concentrations[i])(makeElemArg(_current_elem), 0u));
      nuclide_val = 0.0;
    }
  }

  // We take the max with zero to remove non-physical negatives generated by the mass transport
  // solve in some cases.
  return std::max(val, 0.0);
}
