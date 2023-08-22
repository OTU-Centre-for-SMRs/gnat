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
  const std::vector<VariableName> & _exclude;
};
