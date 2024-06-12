#pragma once

#include "Kernel.h"

// A class which computes a fission heating source term for thermal condution coupled with
// neutronics.
// TODO: Add in the fission heating cross-section to the transport materials section.
class FissionHeatSource : public Kernel
{
public:
  static InputParameters validParams();

  FissionHeatSource(const InputParameters & parameters);

protected:
  virtual Real computeQpResidual() override;

  const unsigned int _num_groups;

  // A vector containing the values of the scalar fluxes, indexed by energy group.
  std::vector<const VariableValue *> _group_scalar_fluxes;

  // A vector containing the fission heating cross-section. Roughly equivalent to \kappa_{g} *
  // \Sigma_{f,g}.
  const ADMaterialProperty<std::vector<Real>> & _kappa_fission;

  // A factor to scale the fission power by.
  const Real _scaling_factor;
}; // class FissionHeatSource
