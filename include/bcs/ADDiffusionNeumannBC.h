#pragma once

#include "ADIntegratedBC.h"

class ADDiffusionNeumannBC : public ADIntegratedBC
{
public:
  static InputParameters validParams();

  ADDiffusionNeumannBC(const InputParameters & parameters);

protected:
  virtual ADReal computeQpResidual() override;

  // Incoming net current for the current group.
  Real _j_g_net;
}; // class ADDiffusionNeumannBC
