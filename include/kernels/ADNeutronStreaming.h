#pragma once

#include "ADNeutronBaseKernel.h"

class ADNeutronStreaming : public ADNeutronBaseKernel
{
public:
  static InputParameters validParams();

  ADNeutronStreaming(const InputParameters & parameters);

protected:
  virtual ADReal computeQpResidual() override;

  const unsigned int _ordinate_index; // n
}; // class ADNeutronStreaming
