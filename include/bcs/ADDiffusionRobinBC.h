#pragma once

#include "ADIntegratedBC.h"

class ADDiffusionRobinBC : public ADIntegratedBC
{
public:
  static InputParameters validParams();

  ADDiffusionRobinBC(const InputParameters & parameters);

protected:
  virtual ADReal computeQpResidual() override;

  // Boundary condition transport correction factor for the current group.
  Real _e_g;
  // Incoming partial current for the current group.
  Real _j_g_inc;
}; // class ADDiffusionRobin
