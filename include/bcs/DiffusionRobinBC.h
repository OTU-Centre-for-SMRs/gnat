#pragma once

#include "IntegratedBC.h"

class DiffusionRobinBC : public IntegratedBC
{
public:
  static InputParameters validParams();

  DiffusionRobinBC(const InputParameters & parameters);

protected:
  virtual Real computeQpResidual() override;
  virtual Real computeQpJacobian() override;

  // Boundary condition transport correction factor for the current group.
  Real _e_g;
  // Incoming partial current for the current group.
  Real _j_g_inc;
}; // class DiffusionRobinBC
