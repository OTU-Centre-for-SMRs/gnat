#pragma once

#include "SNBaseBC.h"

class SNReflectiveBC : public SNBaseBC
{
public:
  static InputParameters validParams();

  SNReflectiveBC(const InputParameters & parameters);

protected:
  static bool vecEquals(const RealVectorValue & first,
                        const RealVectorValue & second,
                        const Real & tol = libMesh::TOLERANCE);

  virtual Real computeQpResidual() override;
  virtual Real computeQpJacobian() override;
  virtual Real computeQpOffDiagJacobian(unsigned int jvar) override;

  unsigned int getReflectedIndex(const RealVectorValue & normal);

  const unsigned int _ordinate_index; // n

  std::set<unsigned int> _jvar_map;
  std::vector<const VariableValue *> _reflected_ordinates;

  const std::vector<RealVectorValue> & _unique_normals;
};
