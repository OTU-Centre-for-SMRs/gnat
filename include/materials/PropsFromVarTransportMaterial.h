#pragma once

#include "EmptyTransportMaterial.h"

class PropsFromVarTransportMaterial : public EmptyTransportMaterial
{
public:
  static InputParameters validParams();

  PropsFromVarTransportMaterial(const InputParameters & parameters);

protected:
  virtual void computeQpProperties() override;

  // Total cross section variables for each group.
  std::vector<const VariableValue *> _sigma_t_g;

  // Scattering matrix variables for each group.
  std::vector<const VariableValue *> _sigma_s_g_prime_g_l;

  // The scattering anisotropy and maximum number of scattering moments.
  const unsigned int _anisotropy;
  const unsigned int _max_moments;

  // Fission neutron production cross sections for each group.
  std::vector<const VariableValue *> _nu_sigma_f_g;

  // Fission neutron spectra values for each group.
  std::vector<const VariableValue *> _chi_f_g;

  // Fission heating values for each group.
  std::vector<const VariableValue *> _heating_g;

  // Particle inverse velocity values for each group.
  std::vector<const VariableValue *> _inv_v_g;

  // Particle diffusion coefficient for each group.
  std::vector<const VariableValue *> _diffusion_g;

  // Absorption cross section for each group.
  std::vector<const VariableValue *> _sigma_a_g;
};
