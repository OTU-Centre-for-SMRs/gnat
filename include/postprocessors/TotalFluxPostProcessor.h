#pragma once

#include "ElementIntegralPostprocessor.h"

// A class that computes the volume and energy integrated scalar flux.
class TotalFluxPostProcessor : public ElementIntegralPostprocessor
{
public:
  static InputParameters validParams();

  TotalFluxPostProcessor(const InputParameters & parameters);

protected:
  virtual Real computeQpIntegral() override;

  // Total number of spectral energy groups.
  const unsigned int _num_groups;

  // The required scalar fluxes.
  std::vector<const VariableValue *> _group_scalar_fluxes;
}; // class TotalFluxPostProcessor
