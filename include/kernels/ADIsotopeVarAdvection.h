#pragma once

#include "ADKernel.h"

class ADIsotopeVarAdvection : public ADKernel
{
public:
  static InputParameters validParams();

  ADIsotopeVarAdvection(const InputParameters & parameters);

protected:
  virtual ADReal computeQpResidual() override;

  // Dummy values in case the problem is a lower dimensionality than 3.
  ADVariableValue _u_vel_temp;
  ADVariableValue _v_vel_temp;
  ADVariableValue _w_vel_temp;

  ADVectorVariableValue _vel_temp;

  // Components of the coupled velocity.
  const ADVariableValue & _u_vel;
  const ADVariableValue & _v_vel;
  const ADVariableValue & _w_vel;

  // Coupled velocity for a VectorValue instead of individual components.
  const ADVectorVariableValue & _vel;

  unsigned int _mesh_dims;
};
