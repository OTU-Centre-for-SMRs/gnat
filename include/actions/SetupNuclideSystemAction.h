#pragma once

#include "Action.h"

class SetupNuclideSystemAction : public Action
{
public:
  static InputParameters validParams();

  SetupNuclideSystemAction(const InputParameters & parameters);

  virtual void act() override;
};
