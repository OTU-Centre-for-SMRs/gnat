#pragma once

#include "SAAFBaseKernel.h"

class SAAFStreaming : public SAAFBaseKernel
{
public:
  static InputParameters validParams();

  SAAFStreaming(const InputParameters & parameters);

protected:
  virtual Real computeQpResidual() override;
  virtual Real computeQpJacobian() override;
}; // class ADSAAFStreaming
