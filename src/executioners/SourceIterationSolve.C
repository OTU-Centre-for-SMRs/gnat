#include "SourceIterationSolve.h"

#include "FEProblem.h"
#include "NonlinearSystem.h"
#include "AuxiliarySystem.h"

InputParameters
SourceIterationSolve::validParams()
{
  InputParameters params = FEProblemSolve::validParams();
  params.addClassDescription(
      "A custom executioner to solve the multi-group transport equation using either: source "
      "iteration or direct solves (PJFNK) for the within-group equations; and Gauss-Seidel or "
      "forward substitution for the group-to-group equations.");
  params.addRequiredRangeCheckedParam<unsigned int>(
      "num_groups", "0<num_groups", "The number of spectral energy groups to use.");
  params.addRequiredRangeCheckedParam<unsigned int>(
      "num_ordinates", "0<num_ordinates", "The number of discrete ordinates to use.");
  params.addRequiredParam<std::vector<NonlinearSystemName>>(
      "nonlinear_system_names",
      "The names of the nonlinear systems for each angular flux variable. For monolithic "
      "within-group solves, the user must provide a nonlinear system for each energy group. For "
      "source-iteration solves, the user must provide a nonlinear system for each energy group and "
      "each flux ordinate. The nonlinear systems are expected to be provided in order of group "
      "first, flux ordinate second.");
  params.addParam<bool>("enable_upscattering",
                        true,
                        "Whether upscattering should be included (using Gauss-Seidel iteration). "
                        "If this parameter is set to false, thermal upscattering is ignored and "
                        "forward substitution is used.");
  params.addParam<MooseEnum>("within_group_solve_type",
                             MooseEnum("monolithic source_iteration", "monolithic"),
                             "The solution strategy to use for the within-group equations.");
  params.addRangeCheckedParam<unsigned int>(
      "max_inner_iterations",
      1000,
      "0<max_inner_iterations",
      "The maximum number of iterations for the within-group equations. Only valid when "
      "source-iteration is selected as the solution scheme.");
  params.addRangeCheckedParam<unsigned int>(
      "max_outer_iterations",
      1000,
      "0<max_outer_iterations",
      "The maximum number of iterations for the group-to-group equations. Only valid when "
      "upscattering is enabled.");
  params.addRangeCheckedParam<Real>(
      "inner_absolute_tolerance",
      1e-8,
      "0.0<inner_absolute_tolerance",
      "The maximum absolute tolerance on the residual of the within-group equations. Only valid "
      "when source-iteration is selected as the solution scheme.");
  params.addRangeCheckedParam<Real>(
      "outer_absolute_tolerance",
      1e-8,
      "0.0<outer_absolute_tolerance",
      "The maximum absolute tolerance on the residual of the group-to-group equations. Only valid "
      "when upscattering is enabled.");

  return params;
}

SourceIterationSolve::SourceIterationSolve(Executioner & ex)
  : FEProblemSolve(ex),
    _num_groups(getParam<unsigned int>("num_groups")),
    _num_ordinates(getParam<unsigned int>("num_ordinates")),
    _enable_thermal_upscattering(getParam<bool>("enable_upscattering")),
    _enable_source_iteration(getParam<MooseEnum>("within_group_solve_type") == "source_iteration"),
    _max_num_inner_iterations(getParam<unsigned int>("max_inner_iterations")),
    _max_num_outer_iterations(getParam<unsigned int>("max_outer_iterations")),
    _inner_absolute_tolerance(getParam<Real>("inner_absolute_tolerance")),
    _outer_absolute_tolerance(getParam<Real>("outer_absolute_tolerance"))
{
  // Fetch the nonlinear systems.
  const auto & sys_names = getParam<std::vector<NonlinearSystemName>>("nonlinear_system_names");
  if (sys_names.size() != _num_groups * _num_ordinates && _enable_source_iteration)
    mooseError("Mismatch between the number of provided nonlinear systems and the number of groups "
               "* number of ordinates.");
  else if (sys_names.size() != _num_groups && !_enable_source_iteration)
    mooseError(
        "Mismatch between the number of provided nonlinear systems and the number of groups.");

  if (_enable_source_iteration)
  {
    _angular_flux_nonlinear_system_numbers.reserve(_num_groups * _num_ordinates);
    _angular_flux_nonlinear_systems.reserve(_num_groups * _num_ordinates);
    for (unsigned int g = 0u; g < _num_groups; ++g)
    {
      for (unsigned int n = 0u; n < _num_ordinates; ++n)
      {
        _angular_flux_nonlinear_system_numbers.emplace_back(
            _problem.nlSysNum(sys_names[g * _num_ordinates + n]));
        _angular_flux_nonlinear_systems.emplace_back(
            &_problem.getNonlinearSystemBase(_angular_flux_nonlinear_system_numbers.back()));
      }
    }
  }
  else
  {
    _angular_flux_nonlinear_system_numbers.reserve(_num_groups);
    _angular_flux_nonlinear_systems.reserve(_num_groups);
    for (unsigned int g = 0u; g < _num_groups; ++g)
    {
      _angular_flux_nonlinear_system_numbers.emplace_back(_problem.nlSysNum(sys_names[g]));
      _angular_flux_nonlinear_systems.emplace_back(
          &_problem.getNonlinearSystemBase(_angular_flux_nonlinear_system_numbers.back()));
    }
  }
}

bool
SourceIterationSolve::solve()
{
  bool converged = false;
  for (MooseIndex(_num_grid_steps) grid_step = 0; grid_step <= _num_grid_steps; ++grid_step)
  {
    if (_enable_thermal_upscattering && _enable_source_iteration)
      converged = gaussSeidelSourceIteration();
    if (_enable_thermal_upscattering && !_enable_source_iteration)
      converged = gaussSeidelFullSolve();
    if (!_enable_thermal_upscattering && _enable_source_iteration)
      converged = forwardSubSourceIteration();
    if (!_enable_thermal_upscattering && !_enable_source_iteration)
      converged = forwardSubFullSolve();

    if (_problem.shouldSolve())
    {
      if (converged)
        _console << COLOR_GREEN << " Solve Converged!" << COLOR_DEFAULT << std::endl;
      else
      {
        _console << COLOR_RED << " Solve Did NOT Converge!" << COLOR_DEFAULT << std::endl;
        return false;
      }
    }
    else
      _console << COLOR_GREEN << " Solve Skipped!" << COLOR_DEFAULT << std::endl;

    if (grid_step != _num_grid_steps)
      _problem.uniformRefine();
  }

  return converged;
}

bool
SourceIterationSolve::gaussSeidelSourceIteration()
{
  mooseError("Gauss-Seidel iteration for the multi-group equations has not been implemented yet.");

  unsigned int num_outer_iterations = 0u;
  Real outer_residual = 1.0;
  // Gauss-Seidel.
  while (_outer_absolute_tolerance < outer_residual &&
         num_outer_iterations < _max_num_outer_iterations)
  {

    num_outer_iterations++;
  }

  return false;
}

bool
SourceIterationSolve::gaussSeidelFullSolve()
{
  mooseError("Gauss-Seidel iteration for the multi-group equations has not been implemented yet.");

  unsigned int num_outer_iterations = 0u;
  Real outer_residual = 1.0;
  // Gauss-Seidel.
  while (_outer_absolute_tolerance < outer_residual &&
         num_outer_iterations < _max_num_outer_iterations)
  {
    num_outer_iterations++;
  }

  return false;
}

bool
SourceIterationSolve::forwardSubSourceIteration()
{
  // Forward substitution.
  unsigned int g = 0u;
  bool should_continue = true;
  while (g < _num_groups && should_continue)
  {
    should_continue = innerIteration(g);
    g++;
  }

  if (!should_continue)
  {
    // Forward substitution failed to converge.
    _console << COLOR_RED << "Forward substitution for group " << g << " did NOT converge!"
             << COLOR_DEFAULT << std::endl;
    return false;
  }
  else
  {
    // Forward substitution converged.
    _console << COLOR_GREEN << "Forward substitution for group " << g << " converged!"
             << COLOR_DEFAULT << std::endl;
    return true;
  }
}

bool
SourceIterationSolve::forwardSubFullSolve()
{
  // Forward substitution.
  unsigned int g = 0u;
  bool should_continue = true;
  while (g < _num_groups && should_continue)
  {
    _angular_flux_nonlinear_systems[g]->residualSetup();
    _problem.solve(_angular_flux_nonlinear_system_numbers[g]);
    should_continue = _problem.nlConverged(_angular_flux_nonlinear_system_numbers[g]);

    g++;
  }

  if (!should_continue)
  {
    // Forward substitution failed to converge.
    _console << COLOR_RED << "Forward substitution for group " << g << " did NOT converge!"
             << COLOR_DEFAULT << std::endl;
    return false;
  }
  else
  {
    // Forward substitution converged.
    _console << COLOR_GREEN << "Forward substitution for group " << g << " converged!"
             << COLOR_DEFAULT << std::endl;
    return true;
  }
}

// The angular fluxes and flux moments MUST be initialized to zero for this to work.
bool
SourceIterationSolve::innerIteration(unsigned int g)
{
  // NonlinearImplicitSystem ExplicitSystem
  // ExplicitSystem sys;

  // TODO: Initial iteration (uncollided flux) required before running the source iteration loop.

  unsigned int num_inner_iterations = 0u;
  Real inner_residual = 1.0;
  bool should_continue = true;
  while (_inner_absolute_tolerance < inner_residual &&
         num_inner_iterations < _max_num_inner_iterations)
  {
    // Check for the stationary point.
    for (unsigned int n = 0u; n < _num_ordinates; ++n)
    {
      inner_residual = std::max(inner_residual,
                                _angular_flux_nonlinear_systems[g * _num_ordinates + n]
                                    ->_initial_residual_after_preset_bcs);
    }

    _console << "Scattering source iteration " << num_inner_iterations << " - Group " << g
             << " angular flux residual maximum norm:\n"
             << COLOR_GREEN << inner_residual << COLOR_DEFAULT << std::endl;

    // Prep and solve for the current iteration.
    for (unsigned int n = 0u; n < _num_ordinates; ++n)
      _angular_flux_nonlinear_systems[g * _num_ordinates + n]->residualSetup();
    for (unsigned int n = 0u; n < _num_ordinates; ++n)
      _problem.solve(_angular_flux_nonlinear_system_numbers[g * _num_ordinates + n]);

    // Check to make sure each flux ordinate converged.
    for (unsigned int n = 0u; n < _num_ordinates; ++n)
      should_continue =
          should_continue &&
          _problem.nlConverged(_angular_flux_nonlinear_system_numbers[g * _num_ordinates + n]);

    // Quit if one of the nonlinear systems fails to converge.
    if (!should_continue)
      break;

    // Flux moments are computed on EXEC_LINEAR and EXEC_TIMESTEP_END by default. In the case of
    // source iteration we need to control when the flux moments are computed so we can
    // appropriately lag the scattering source. Hence, we assume that the exec flag for the flux
    // moments is set to EXEC_CUSTOM and call computeAuxiliaryKernels with that flag.
    _problem.computeAuxiliaryKernels(EXEC_CUSTOM);

    num_inner_iterations++;
  }

  if (inner_residual >= _inner_absolute_tolerance || !should_continue)
  {
    // Source iteration failed to converge.
    _console << COLOR_RED << "Scattering source iteration for group " << g << " did NOT converge!"
             << COLOR_DEFAULT << std::endl;
    return false;
  }
  else
  {
    // Source iteration converged.
    _console << COLOR_GREEN << "Scattering source iteration for group " << g << " converged!"
             << COLOR_DEFAULT << std::endl;
    return true;
  }
}
