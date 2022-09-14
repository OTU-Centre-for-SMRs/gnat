#pragma once

#include "Action.h"

class CommonGnatAction : public Action
{
public:
  static InputParameters validParams();

  CommonGnatAction(const InputParameters & params);

  virtual void act() override;
}; // class CommonGnatAction
