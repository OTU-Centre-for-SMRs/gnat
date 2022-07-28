#pragma once

#include "SNDiracKernelBase.h"

class SAAFDiracKernelBase : public SNDiracKernelBase
{
public:
  static InputParameters validParams();

  SAAFDiracKernelBase(const InputParameters & parameters);

protected:
  Real maxVertexSeparation();
  Real computeQPTau();
  Real computeQPTests();

  const unsigned int _ordinate_index; // n
  const unsigned int _group_index; // g

  const ADMaterialProperty<std::vector<Real>> & _sigma_r_g;
}; // class SAAFDiracKernelBase
