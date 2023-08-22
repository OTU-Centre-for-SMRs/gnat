#pragma once

#include "MooseObjectAction.h"

#include "GnatBase.h"

class AddFVNuclideBCAction : public MooseObjectAction
{
public:
  static InputParameters validParams();

  AddFVNuclideBCAction(const InputParameters & parameters);

  virtual void act() override;

protected:
};
