#pragma once

#include "EmptyTransportMaterial.h"

class AbsorbingTransportMaterial : public EmptyTransportMaterial
{
public:
  static InputParameters validParams();

  AbsorbingTransportMaterial(const InputParameters & parameters);

protected:
  virtual void computeQpProperties() override;

  std::vector<Real> _inv_v_g;
  std::vector<Real> _sigma_t_g;
  std::vector<Real> _diffusion_g;
}; // class AbsorbingTransportMaterial
