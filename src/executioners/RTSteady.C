#include "RTSteady.h"

registerMooseObject("GnatApp", RTSteady);

InputParameters
RTSteady::validParams()
{
  InputParameters params = Steady::validParams();
  params += SourceIterationSolve::validParams();
  params.addClassDescription(
      "A custom executioner to solve the steady-state multi-group transport equation using either: "
      "source iteration or direct solves (PJFNK) for the within-group equations; and Gauss-Seidel "
      "or forward substitution for the group-to-group equations. This overrides some functionality "
      "in the default 'Steady' executioner provided by MOOSE to enabled these iterative radiation "
      "transport schemes.");

  return params;
}

RTSteady::RTSteady(const InputParameters & parameters)
  : Steady(parameters), _source_iteration(*this)
{
  _fixed_point_solve->setInnerSolve(_source_iteration);
}
