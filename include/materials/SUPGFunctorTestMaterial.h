#pragma once

#include "FunctorMaterial.h"

// A material to provide the required properties for the nuclide kernels during testing.
class SUPGFunctorTestMaterial : public FunctorMaterial
{
public:
  static InputParameters validParams();

  SUPGFunctorTestMaterial(const InputParameters & parameters);

protected:
  const std::string & _var_name;

  const Real _diff;
  const RealVectorValue _vel;
}; // class SUPGFunctorTestMaterial
