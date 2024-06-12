#pragma once

#include "ElementIntegralPostprocessor.h"

// A class that computes the total fission power.
class FissionPowerPostProcessor : public ElementIntegralPostprocessor
{
public:
  static InputParameters validParams();

  FissionPowerPostProcessor(const InputParameters & parameters);

protected:
  virtual Real computeQpIntegral() override;

  // Total number of spectral energy groups.
  const unsigned int _num_groups;

  // The required scalar fluxes.
  std::vector<const VariableValue *> _group_scalar_fluxes;

  // A vector containing the fission heating cross-section. Roughly equivalent to \kappa_{g} *
  // \Sigma_{f,g}.
  const ADMaterialProperty<std::vector<Real>> & _kappa_fission;
}; // class FissionPowerPostProcessor
