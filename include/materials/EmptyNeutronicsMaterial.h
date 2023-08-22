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

  // Speed of light in cm s^{-1}.
  static constexpr Real _c_cm = 2.99792458e10;
  static constexpr Real _inv_c_cm = 3.335640952e-11;
  // Speed of light in m s^{-1}.
  static constexpr Real _c_m = 2.997924580e8;
  static constexpr Real _inv_c_m = 3.335640952e-9;

  // A scale factor so we can convert from cm^{-1} to x^{-1} and cm s^{-1} to x s^{-1}
  // (where x is an arbitrary length unit).
  const Real _scaling;

  const unsigned int _num_groups;
  const Particletype _particle;

  const bool _is_saaf;
  const bool _is_diffusion;
  const bool _has_fission;

  // Material properties that neutronics materials are expected to provide.
  ADMaterialProperty<std::vector<Real>> & _mat_inv_v_g;
  ADMaterialProperty<std::vector<Real>> & _mat_sigma_t_g;
  ADMaterialProperty<std::vector<Real>> & _mat_surface_source;
  ADMaterialProperty<std::vector<Real>> & _mat_sigma_s_g_prime_g_l;
  MaterialProperty<unsigned int> & _mat_anisotropy;
  ADMaterialProperty<std::vector<Real>> & _mat_source_moments;
  MaterialProperty<unsigned int> & _mat_src_anisotropy;

  // Material properties for diffusion schemes.
  ADMaterialProperty<std::vector<Real>> * _mat_sigma_r_g;
  ADMaterialProperty<std::vector<Real>> * _mat_diffusion_g;

  // Material properties for fission.
  ADMaterialProperty<std::vector<Real>> * _mat_nu_sigma_f_g;
  ADMaterialProperty<std::vector<Real>> * _mat_chi_f_g;

  // SAAF stabilization parameters.
  ADMaterialProperty<std::vector<Real>> * _mat_saaf_tau;
  Real _saaf_eta;
  Real _saaf_c;
}; // class EmptyNeutronicsMaterial
