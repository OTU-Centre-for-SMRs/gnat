#pragma once

#include "MooseObjectAction.h"

#include "GnatBase.h"

class AddTransportMaterialAction : public MooseObjectAction
{
public:
  static InputParameters validParams();

  AddTransportMaterialAction(const InputParameters & parameters);

  virtual void act() override;

protected:
  std::string _parent_transport_system;
  unsigned int _num_groups;
  MooseEnum _particle;
  MooseEnum _scheme;
  bool _disable_fission;

  bool _is_init;
}; // class AddTransportMaterialAction
