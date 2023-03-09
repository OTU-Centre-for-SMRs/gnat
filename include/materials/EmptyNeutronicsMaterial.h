#pragma once

#include "Material.h"

#include "GnatBase.h"

class EmptyNeutronicsMaterial : public Material
{
public:
  static InputParameters validParams();

  EmptyNeutronicsMaterial(const InputParameters & parameters);

protected:
  virtual void computeQpProperties() override;

  const unsigned int _num_groups;
  const Particletype _particle;

  // Material properties that neutronics materials are expected to provide.
  ADMaterialProperty<std::vector<Real>> & _mat_inv_v_g;
  ADMaterialProperty<std::vector<Real>> & _mat_sigma_t_g;
  ADMaterialProperty<std::vector<Real>> & _mat_sigma_r_g;
  ADMaterialProperty<std::vector<Real>> & _mat_diffusion_g;
  ADMaterialProperty<std::vector<Real>> & _mat_surface_source;
  ADMaterialProperty<std::vector<Real>> & _mat_sigma_s_g_prime_g_l;
  MaterialProperty<unsigned int> & _mat_anisotropy;
  ADMaterialProperty<std::vector<Real>> & _mat_source_moments;
  MaterialProperty<unsigned int> & _mat_src_anisotropy;

  // SAAF stabilization parameters.
  ADMaterialProperty<Real> & _mat_saaf_eta;
  ADMaterialProperty<Real> & _mat_saaf_c;
  Real _saaf_eta;
  Real _saaf_c;
}; // class EmptyNeutronicsMaterial
