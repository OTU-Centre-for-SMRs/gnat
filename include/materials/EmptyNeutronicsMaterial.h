#pragma once

#include "Material.h"

class EmptyNeutronicsMaterial : public Material
{
public:
  static InputParameters validParams();

  EmptyNeutronicsMaterial(const InputParameters & parameters);

protected:
  virtual void computeQpProperties() override;

  const unsigned int _num_groups;

  // Material properties that neutronics materials are expected to provide.
  ADMaterialProperty<std::vector<Real>> & _mat_v_g;
  ADMaterialProperty<std::vector<Real>> & _mat_sigma_r_g;
  ADMaterialProperty<std::vector<Real>> & _mat_surface_source;
  ADMaterialProperty<std::vector<Real>> & _mat_sigma_s_g_prime_g_l;
  MaterialProperty<unsigned int> & _mat_anisotropy;
  ADMaterialProperty<std::vector<Real>> & _mat_source_moments;
  MaterialProperty<unsigned int> & _mat_src_anisotropy;
}; // class EmptyNeutronicsMaterial
