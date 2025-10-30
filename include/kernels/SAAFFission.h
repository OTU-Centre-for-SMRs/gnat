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

  Real computeScalarFlux(unsigned int g_prime);

  // Total number of spectral energy groups.
  const unsigned int _num_groups;

  /*
   * We assume that the vector of flux ordinates is stored in order of group
   * first, direction second. An example for 2 energy groups (G = 2) and a
   * quadrature set with 2 elements (N = 2) is given below:
   *
   * The flux ordinates are indexed as Psi_{g, n}:
   * _group_flux_ordinates[0] = Psi_{1, 1}
   * _group_flux_ordinates[0] = Psi_{1, 2}
   * _group_flux_ordinates[0] = Psi_{2, 1}
   * _group_flux_ordinates[0] = Psi_{2, 2}
   */
  std::map<unsigned int, std::pair<unsigned int, unsigned int>> _jvar_map;
  std::vector<const VariableValue *> _group_flux_ordinates;

  // The neutron production cross-sections.
  const ADMaterialProperty<std::vector<Real>> & _nu_sigma_f_g;
  // The fission production spectra.
  const ADMaterialProperty<std::vector<Real>> & _chi_g;
}; // class SAAFFission
