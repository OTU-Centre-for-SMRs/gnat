#pragma once

#include "FEProblemSolve.h"

class NonlinearSystemBase;

// A class which implements both source iteration for the within-group equations and Gauss-Seidel
// iteration for the multi-group problem with thermal particle upscattering.
class SourceIterationSolve : public FEProblemSolve
{
public:
  static InputParameters validParams();

  SourceIterationSolve(Executioner & ex);

  virtual bool solve() override;

  virtual void setInnerSolve(SolveObject &) override
  {
    mooseError("Cannot set inner solve for SourceIterationSolve.");
  }

protected:
  // The four functions below return true if the solve converges.
  // Solve the multi-group problem with Gauss-Seidel iteration to handle thermal upscattering.
  // Source iteration is used to resolve the within-group equations.
  bool gaussSeidelSourceIteration();
  // Solve the multi-group problem with Gauss-Seidel iteration to handle thermal upscattering.
  // A full matrix solve is used to resolve the within-group equations.
  bool gaussSeidelFullSolve();
  // Solve the multi-group problem using forward substitution, making the assumption that thermal
  // upscattering can be ignored. Source iteration is used to resolve the within-group equations.
  bool forwardSubSourceIteration();
  // Solve the multi-group problem using forward substitution, making the assumption that thermal
  // upscattering can be ignored. A full matrix solve is used to resolve the within-group equations.
  bool forwardSubFullSolve();
  // The source iteration solver.
  bool innerIteration(unsigned int g);

  std::vector<unsigned int> _angular_flux_nonlinear_system_numbers;
  std::vector<NonlinearSystemBase *> _angular_flux_nonlinear_systems;

  // Number of energy groups and number of discrete ordinates, respectively.
  const unsigned int _num_groups;
  const unsigned int _num_ordinates;

  // Whether thermal upscattering should be enabled for the multi-group problem or not.
  // If not enabled, forward substitution (starting with the high-energy groups) is used
  // to resolve the multigroup problem.
  const bool _enable_thermal_upscattering;

  // Whether source iteration should be used or not. If source iteration is not enabled the
  // within-group equations are solved together.
  const bool _enable_source_iteration;

  // Maximum number of source iterations.
  const unsigned int _max_num_inner_iterations;
  // Maximum number of Gauss-Seidel (GS) iterations for upscattering multigroup problems.
  const unsigned int _max_num_outer_iterations;

  // Source iteration convergence criteria.
  const Real _inner_absolute_tolerance;
  // Gauss-Seidel convergence criteria.
  const Real _outer_absolute_tolerance;
}; // class SourceIterationSolve
