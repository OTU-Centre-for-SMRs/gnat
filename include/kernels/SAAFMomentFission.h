#pragma once

#include "SAAFBaseKernel.h"

// A class to compute the residual contribution from fission for the neutron-specific form of the
// discrete ordinates radiation transport equation. This kernel operates on scalar flux moments and
// does not assemble the off-diagonal Jacobian.
class SAAFMomentFission : public SAAFBaseKernel
{
public:
  static InputParameters validParams();

  SAAFMomentFission(const InputParameters & parameters);

protected:
  virtual Real computeQpResidual() override;
  virtual Real computeQpJacobian() override;

  // Total number of spectral energy groups.
  const unsigned int _num_groups;

  // The required scalar fluxes.
  std::vector<const VariableValue *> _group_scalar_fluxes;

  // The neutron production cross-sections.
  const ADMaterialProperty<std::vector<Real>> & _nu_sigma_f_g;
  // The fission production spectra.
  const ADMaterialProperty<std::vector<Real>> & _chi_g;
}; // class SAAFMomentFission
