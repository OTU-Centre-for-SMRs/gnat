#pragma once

#include "Kernel.h"

class DiffusionFission : public Kernel
{
public:
  static InputParameters validParams();

  DiffusionFission(const InputParameters & parameters);

protected:
  virtual Real computeQpResidual() override;
  virtual Real computeQpJacobian() override;
  virtual Real computeQpOffDiagJacobian(unsigned int jvar) override;

  // The current group (g) and the number of spectral energy groups (G).
  const unsigned int _group_index;
  const unsigned int _num_groups;

  // The required scalar fluxes.
  std::map<unsigned int, unsigned int> _jvar_map;
  std::vector<const VariableValue *> _group_scalar_fluxes;

  // The neutron production cross-sections.
  const ADMaterialProperty<std::vector<Real>> & _nu_sigma_f_g;
  // The fission production spectra.
  const ADMaterialProperty<std::vector<Real>> & _chi_g;
}; // class ADDiffusionFission
