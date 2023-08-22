#pragma once

#include "SNBaseBC.h"

// A boundary condition which imposes a surface current along the neutron direction of travel
// aligned with the surface normal. Equivalent to the angular distribution being a Dirac delta
// function shifted by the surface normal: $\delta(\hat{\Omega} - \hat{n})$.
class SNNormalCurrentBC : public SNBaseBC
{
public:
  static InputParameters validParams();

  SNNormalCurrentBC(const InputParameters & parameters);

protected:
  virtual Real computeQpResidual() override;
  virtual Real computeQpJacobian() override;

  const unsigned int _num_groups;     // G
  const unsigned int _group_index;    // g
  const unsigned int _ordinate_index; // n

  // The incoming current.
  Real _current;

  // Anisotropy of the imposed current. Used to increase the accuracy of the SH expansion.
  const unsigned int _current_anisotropy;
}; // class ADSNNormalCurrentBC
