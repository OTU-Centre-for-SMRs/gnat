#pragma once

#include "EmptyNeutronicsMaterial.h"

class VoidNeutronicsMaterial : public EmptyNeutronicsMaterial
{
public:
  static InputParameters validParams();

  VoidNeutronicsMaterial(const InputParameters & parameters);

protected:
  virtual void computeQpProperties() override;

  std::vector<Real> _v_g;
}; // VoidNeutronicsMaterial
