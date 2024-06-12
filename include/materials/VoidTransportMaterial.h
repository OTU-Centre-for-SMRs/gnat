#pragma once

#include "EmptyTransportMaterial.h"

class VoidTransportMaterial : public EmptyTransportMaterial
{
public:
  static InputParameters validParams();

  VoidTransportMaterial(const InputParameters & parameters);

protected:
  virtual void computeQpProperties() override;

  std::vector<Real> _inv_v_g;
  std::vector<Real> _diffusion_g;
}; // VoidTransportMaterial
