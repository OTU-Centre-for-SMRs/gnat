#pragma once

#include "Transient.h"
#include "SourceIterationSolve.h"

// An executioner which overrides some functionality in the default transient MOOSE executioner to
// enable iterative solutions to the radiation transport equation.
class RTTransient : public Transient
{
public:
  static InputParameters validParams();

  RTTransient(const InputParameters & parameters);

protected:
  SourceIterationSolve _source_iteration;
}; // class RTTransient
