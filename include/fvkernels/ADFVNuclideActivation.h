#pragma once

#include "FVElementalKernel.h"

class ADFVNuclideActivation : public FVElementalKernel
{
public:
  static InputParameters validParams();

  ADFVNuclideActivation(const InputParameters & parameters);

protected:
  virtual ADReal computeQpResidual() override;

  // Number of energy groups.
  const unsigned int _num_groups;

  // A vector to store all of the group scalar fluxes.
  std::vector<const ADVariableValue *> _group_scalar_fluxes;

  // A vector to store all of the isotope number densities which form the
  // current isotope under neutron bombardment.
  std::vector<const Moose::Functor<ADReal> *> _isotope_number_densities;

  // A vector to store all of the group activation cross-sections.
  std::vector<Real> _sigma_a_g;
}; // class ADFVNuclideActivation
