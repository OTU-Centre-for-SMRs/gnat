#pragma once

#include "IntegratedBC.h"

class DiffusionNeumannBC : public IntegratedBC
{
public:
  static InputParameters validParams();

  DiffusionNeumannBC(const InputParameters & parameters);

protected:
  virtual Real computeQpResidual() override;

  // Incoming net current for the current group.
  Real _j_g_net;
}; // class ADDiffusionNeumannBC
