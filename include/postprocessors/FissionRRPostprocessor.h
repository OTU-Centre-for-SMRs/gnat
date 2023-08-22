#pragma once

#include "ElementIntegralPostprocessor.h"

// A class that computes the integrated fission reaction rate.
class FissionRRPostprocessor : public ElementIntegralPostprocessor
{
public:
  static InputParameters validParams();

  FissionRRPostprocessor(const InputParameters & parameters);

protected:
  virtual Real computeQpIntegral() override;

  // Total number of spectral energy groups.
  const unsigned int _num_groups;

  // The required scalar fluxes.
  std::vector<const VariableValue *> _group_scalar_fluxes;

  // The neutron production cross-sections.
  const ADMaterialProperty<std::vector<Real>> & _nu_sigma_f_g;
}; // class FissionRRPostprocessor
