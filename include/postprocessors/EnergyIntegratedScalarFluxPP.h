#pragma once

#include "ElementIntegralPostprocessor.h"

// A class that computes the energy integrated scalar particle flux.
class EnergyIntegratedScalarFluxPP : public ElementIntegralPostprocessor
{
public:
  static InputParameters validParams();

  EnergyIntegratedScalarFluxPP(const InputParameters & parameters);

protected:
  virtual Real computeQpIntegral() override;

  // Total number of spectral energy groups.
  const unsigned int _num_groups;

  // The required scalar fluxes.
  std::vector<const VariableValue *> _group_scalar_fluxes;
}; // class EnergyIntegratedScalarFluxPP
