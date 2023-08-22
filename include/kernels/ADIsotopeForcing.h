#pragma once

#include "ADIsotopeBase.h"

class Function;

// A forcing function for testing isotope kernels.
class ADIsotopeForcing : public ADIsotopeBase
{
public:
  static InputParameters validParams();

  ADIsotopeForcing(const InputParameters & parameters);

protected:
  virtual ADReal computeQpResidual() override;

  const Function & _forcing;
}; // class ADIsotopeForcing
