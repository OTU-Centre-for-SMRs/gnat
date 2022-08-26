#pragma once

#include "ADIsotopeBase.h"

class ADIsotopeDepletion : public ADIsotopeBase
{
public:
  static InputParameters validParams();

  ADIsotopeDepletion(const InputParameters & parameters);

protected:
  virtual ADReal computeQpResidual() override;

  // Number of energy groups.
  const unsigned int _num_groups;

  // A vector to store all of the group scalar fluxes.
  std::vector<const ADVariableValue *> _group_scalar_fluxes;

  // A vector to store all of the group absorption cross-sections.
  std::vector<Real> _sigma_a_g;
}; // class ADIsotopeDepletion
