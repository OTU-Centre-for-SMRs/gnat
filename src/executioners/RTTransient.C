#include "RTTransient.h"

registerMooseObject("GnatApp", RTTransient);

InputParameters
RTTransient::validParams()
{
  InputParameters params = Transient::validParams();
  params += SourceIterationSolve::validParams();
  params.addClassDescription(
      "A custom executioner to solve the transient multi-group transport equation using either: "
      "source iteration or direct solves (PJFNK) for the within-group equations; and Gauss-Seidel "
      "or forward substitution for the group-to-group equations. This overrides some functionality "
      "in the default 'Transient' executioner provided by MOOSE to enabled these iterative "
      "radiation transport schemes.");

  return params;
}

RTTransient::RTTransient(const InputParameters & parameters)
  : Transient(parameters), _source_iteration(*this)
{
  _fixed_point_solve->setInnerSolve(_source_iteration);
}
