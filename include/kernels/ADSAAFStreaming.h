#pragma once

#include "ADSAAFBaseKernel.h"

class ADSAAFStreaming : public ADSAAFBaseKernel
{
public:
  static InputParameters validParams();

  ADSAAFStreaming(const InputParameters & parameters);

protected:
  virtual ADReal computeQpResidual() override;
}; // class ADSAAFStreaming
