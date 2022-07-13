#pragma once

#include "Material.h"

class AbsorbingNeutronicsMaterial : public Material
{
public:
  static InputParameters validParams();

  AbsorbingNeutronicsMaterial(const InputParameters & parameters);

protected:
  virtual void computeQpProperties() override;

  const unsigned int _num_groups;

  std::vector<Real> _v_g;
  std::vector<Real> _sigma_r_g;

  ADMaterialProperty<std::vector<Real>> & _mat_v_g;
  ADMaterialProperty<std::vector<Real>> & _mat_sigma_r_g;
}; // class AbsorbingNeutronicsMaterial
