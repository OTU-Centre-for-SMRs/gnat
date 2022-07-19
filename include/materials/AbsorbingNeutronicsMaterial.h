#pragma once

#include "EmptyNeutronicsMaterial.h"

class AbsorbingNeutronicsMaterial : public EmptyNeutronicsMaterial
{
public:
  static InputParameters validParams();

  AbsorbingNeutronicsMaterial(const InputParameters & parameters);

protected:
  virtual void computeQpProperties() override;

  std::vector<Real> _v_g;
  std::vector<Real> _sigma_r_g;
}; // class AbsorbingNeutronicsMaterial
