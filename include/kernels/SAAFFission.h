#pragma once

#include "SAAFBaseKernel.h"

// A class to compute the residual contribution from fission for the neutron-specific form of the
// discrete ordinates radiation transport equation. This kernel assembles the full fission Jacobian.
class SAAFFission : public SAAFBaseKernel
{
public:
  static InputParameters validParams();

  SAAFFission(const InputParameters & parameters);

protected:
  virtual Real computeQpResidual() override;
  virtual Real computeQpJacobian() override;
  virtual Real computeQpOffDiagJacobian(unsigned int jvar) override;
  // Total number of spectral energy groups.
  const unsigned int _num_groups;

  // A map between the coupleable variable ID and the ordinate / group index.
  std::map<unsigned int, std::pair<unsigned int, unsigned int>> _jvar_map;

  // The required scalar fluxes which we actually use for calculations.
  std::vector<const VariableValue *> _group_scalar_fluxes;

  // The neutron production cross-sections.
  const ADMaterialProperty<std::vector<Real>> & _nu_sigma_f_g;
  // The fission production spectra.
  const ADMaterialProperty<std::vector<Real>> & _chi_g;
}; // class SAAFFission
