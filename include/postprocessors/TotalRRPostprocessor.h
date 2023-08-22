#pragma once

#include "ElementIntegralPostprocessor.h"

// A class that computes the integrated total reaction rate.
class TotalRRPostprocessor : public ElementIntegralPostprocessor
{
public:
  static InputParameters validParams();

  TotalRRPostprocessor(const InputParameters & parameters);

protected:
  virtual Real computeQpIntegral() override;

  // Total number of spectral energy groups.
  const unsigned int _num_groups;

  // The required scalar fluxes.
  std::vector<const VariableValue *> _group_scalar_fluxes;

  const ADMaterialProperty<std::vector<Real>> & _sigma_t_g;
}; // class TotalRRPostprocessor
