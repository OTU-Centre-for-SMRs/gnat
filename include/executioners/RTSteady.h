#pragma once

#include "Steady.h"
#include "SourceIterationSolve.h"

// An executioner which overrides some functionality in the default steady-state MOOSE executioner
// to enable iterative solutions to the radiation transport equation.
class RTSteady : public Steady
{
public:
  static InputParameters validParams();

  RTSteady(const InputParameters & parameters);

protected:
  SourceIterationSolve _source_iteration;
}; // class RTSteady
