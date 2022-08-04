#pragma once

#include "SNBaseDiracKernel.h"

class SAAFBaseDiracKernel : public SNBaseDiracKernel
{
public:
  static InputParameters validParams();

  SAAFBaseDiracKernel(const InputParameters & parameters);

protected:
  Real maxVertexSeparation();
  Real computeQPTau();
  Real computeQPTests();

  const unsigned int _ordinate_index; // n
  const unsigned int _group_index; // g

  const ADMaterialProperty<std::vector<Real>> & _sigma_r_g;

  // SAAF stabilization parameters.
  const ADMaterialProperty<Real> & _saaf_eta;
  const ADMaterialProperty<Real> & _saaf_c;
}; // class SAAFBaseDiracKernel
