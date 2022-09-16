#pragma once

#include "Action.h"

class SetupIsotopeSystemAction : public Action
{
public:
  static InputParameters validParams();

  SetupIsotopeSystemAction(const InputParameters & parameters);

  virtual void act() override;
};
