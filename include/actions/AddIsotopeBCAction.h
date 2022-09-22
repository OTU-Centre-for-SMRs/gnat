#pragma once

#include "MooseObjectAction.h"

#include "GnatBase.h"

class AddIsotopeBCAction : public MooseObjectAction
{
public:
  static InputParameters validParams();

  AddIsotopeBCAction(const InputParameters & parameters);

  virtual void act() override;

protected:
  void applyIsotopeParameters(InputParameters & params);

  // The coordinate system type and dimensionality.
  ProblemType _p_type;

  const std::vector<VariableName> & _master_isotope_list;
};
