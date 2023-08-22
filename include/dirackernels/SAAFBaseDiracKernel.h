#pragma once

#include "SNBaseDiracKernel.h"

class SAAFBaseDiracKernel : public SNBaseDiracKernel
{
public:
  static InputParameters validParams();

  SAAFBaseDiracKernel(const InputParameters & parameters);

protected:
  Real computeQpTests();

  const unsigned int _ordinate_index; // n
  const unsigned int _group_index;    // g

  const ADMaterialProperty<std::vector<Real>> & _sigma_r_g;

  // SAAF stabilization parameters.
  const ADMaterialProperty<std::vector<Real>> & _saaf_tau;
}; // class SAAFBaseDiracKernel
