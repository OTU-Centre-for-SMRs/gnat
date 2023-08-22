#pragma once

#include "Kernel.h"

class SASFAdvection : public Kernel
{
public:
  static InputParameters validParams();

  SASFAdvection(const InputParameters & parameters);

protected:
  virtual Real computeQpResidual() override;
  virtual Real computeQpJacobian() override;

  // The source location.
  const Point _source_location;

  // Divergence multiplicity.
  Real _div_mult;

  // g
  const unsigned int _group_index;
  // Total cross-section.
  const ADMaterialProperty<std::vector<Real>> & _sigma_t_g;
  // SAAF stabilization parameters.
  const ADMaterialProperty<std::vector<Real>> & _saaf_tau;
}; // class SASFAdvection
