#pragma once

#include "ADSNBaseKernel.h"

class ADDFEMStreaming : public ADSNBaseKernel
{
public:
  static InputParameters validParams();

  ADDFEMStreaming(const InputParameters & parameters);

protected:
  virtual ADReal computeQpResidual() override;

  const unsigned int _ordinate_index; // n
}; // class ADDFEMStreaming
