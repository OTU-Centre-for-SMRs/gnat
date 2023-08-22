#include "AddFVNuclideBCAction.h"

#include "FEProblem.h"
#include "MobileDepletionSystemAction.h"
#include "DepletionLibraryAction.h"
#include "NSFVAction.h"

#include "Nuclide.h"

registerMooseAction("GnatApp", AddFVNuclideBCAction, "add_fv_bc");

InputParameters
AddFVNuclideBCAction::validParams()
{
  auto params = MooseObjectAction::validParams();

  params.addClassDescription("Adds a FV boundary condition to all isotopes in the isotope system.");
  params.addParam<bool>("include_full_chain",
                        false,
                        "Whether this boundary condition should apply to the entire depletion "
                        "chain, or just the initial elements.");

  return params;
}

AddFVNuclideBCAction::AddFVNuclideBCAction(const InputParameters & parameters)
  : MooseObjectAction(parameters)
{
}

void
AddFVNuclideBCAction::act()
{
  const auto nuclide_actions = _awh.getActions<MobileDepletionSystemAction>();
  if (nuclide_actions.size() != 1u)
    mooseError("There can only be a single MobileDepletionSystemAction at a time. There are currently " +
               Moose::stringify(nuclide_actions.size() + " MobileDepletionSystemActions."));
  const auto & nuclide_action = (*nuclide_actions[0u]);

  // Fetch the required coupled depletion library.
  const DepletionLibraryAction * coupled_depletion_lib = nullptr;
  {
    const auto depletion_lib_actions = _awh.getActions<DepletionLibraryAction>();
    if (depletion_lib_actions.size() != 1u)
      mooseError("The input file must have a single Depletion Library. There are currently " +
                 Moose::stringify(depletion_lib_actions.size() + " DepletionLibraryActions."));
    coupled_depletion_lib = depletion_lib_actions[0u];
  }

  // Fetch the coupled Navier-Stokes finite volume action.
  const NSFVAction * coupled_ns_fv = nullptr;
  if (nuclide_action.getParam<bool>("using_moose_ns_fv"))
  {
    const auto ns_fv_actions = _awh.getActions<NSFVAction>();
    if (ns_fv_actions.size() != 1u)
      mooseError("The input file must have a single Navier-Stokes system. There are currently " +
                 Moose::stringify(ns_fv_actions.size() + " NSFVActions."));
    coupled_ns_fv = ns_fv_actions[0u];
  }

  std::unordered_map<std::string, Real> total_nuclide_list;
  {
    const auto & elements = nuclide_action.getParam<std::vector<std::string>>("elements");
    const auto & element_atom_fractions =
        nuclide_action.getParam<std::vector<Real>>("element_atom_fractions");
    const auto & extra_nuclides =
        nuclide_action.getParam<std::vector<std::string>>("extra_nuclides");
    const auto & extra_nuclide_atom_fractions =
        nuclide_action.getParam<std::vector<Real>>("extra_nuclide_atom_fractions");

    // First store everything in atom fractions.
    for (unsigned int i = 0u; i < elements.size(); ++i)
    {
      auto abundances = NuclearData::Nuclide::getAbundances(elements[i]);
      for (auto & [nuclide, abundance] : abundances)
      {
        if (total_nuclide_list.count(nuclide) == 0u)
          total_nuclide_list.emplace(nuclide, abundance * element_atom_fractions[i]);
        else
          _console << COLOR_YELLOW << "The nuclide " << nuclide
                   << " has been provided multiple times. This iteration (element atom fraction of "
                   << element_atom_fractions[i]
                   << ") will be ignored when computing the composition of the mixture.\n"
                   << COLOR_DEFAULT;
      }
    }

    for (unsigned int i = 0u; i < extra_nuclides.size(); ++i)
    {
      if (total_nuclide_list.count(extra_nuclides[i]) == 0u)
        total_nuclide_list.emplace(extra_nuclides[i], extra_nuclide_atom_fractions[i]);
      else
        _console << COLOR_YELLOW << "The nuclide " << extra_nuclides[i]
                 << " has been provided multiple times. This iteration (nuclide atom fraction of "
                 << extra_nuclide_atom_fractions[i]
                 << ") will be ignored when computing the composition of the mixture.\n"
                 << COLOR_DEFAULT;
    }

    // Now perform a normalization to compute weight fractions.
    {
      Real total_weight = 0.0;
      for (auto & [nuclide, fraction] : total_nuclide_list)
      {
        fraction *= NuclearData::Nuclide::getAtomicMass(nuclide);
        total_weight += fraction;
      }
      for (auto & [nuclide, fraction] : total_nuclide_list)
        fraction /= total_weight;
    }

    // Sanity check to make sure everything adds up to 1.0. If not, warn the user and scale each
    // weight fraction to ensure the results are conservative.
    {
      Real sum = 0.0;
      for (auto & [nuclide, fraction] : total_nuclide_list)
        sum += fraction;

      if (!MooseUtils::absoluteFuzzyEqual(sum, 1.0))
      {
        _console << COLOR_YELLOW << "The sum of all computed weight fractions (" << sum
                 << ") is not 1. Each weight fraction will be multiplied by " << (1.0 / sum)
                 << " to remain conservative.\n"
                 << COLOR_DEFAULT;
        for (auto & [nuclide, fraction] : total_nuclide_list)
          fraction /= sum;
      }
    }
  }

  std::vector<std::string> nuclides;
  nuclides.reserve(total_nuclide_list.size());
  for (const auto & [nuclide, density] : total_nuclide_list)
    nuclides.emplace_back(nuclide);

  // Fetch the entire depletion chain.
  coupled_depletion_lib->getNuclidesFromInitial(nuclides);
  for (const auto & nuclide : nuclides)
  {
    if (total_nuclide_list.count(nuclide) == 0u)
      total_nuclide_list.emplace(nuclide, 0.0);
  }

  // Remove nuclides that do not exist in the depletion system.
  {
    std::vector<std::string> remove;
    for (const auto & [nuclide, weight] : total_nuclide_list)
    {
      if (!coupled_depletion_lib->hasNuclide(nuclide))
        remove.emplace_back(nuclide);
    }

    for (const auto & nuclide : remove)
      total_nuclide_list.erase(nuclide);
  }

  for (const auto & [nuclide, fraction] : total_nuclide_list)
  {
    _moose_object_pars.set<NonlinearVariableName>("variable") = nuclide + "_mass_fraction";
    _moose_object_pars.set<Real>("inflow_rate") = fraction;

    if (coupled_ns_fv)
    {
      _moose_object_pars.set<MooseFunctorName>("density") =
          coupled_ns_fv->getParam<MooseFunctorName>("density");

      _moose_object_pars.set<MooseFunctorName>("u") = "vel_x";
      if (_problem->mesh().dimension() >= 2u)
        _moose_object_pars.set<MooseFunctorName>("v") = "vel_y";
      if (_problem->mesh().dimension() >= 3u)
        _moose_object_pars.set<MooseFunctorName>("w") = "vel_z";
    }
    else
    {
      _moose_object_pars.set<MooseFunctorName>("density") =
          nuclide_action.getParam<MooseFunctorName>("density");

      _moose_object_pars.set<MooseFunctorName>("u") =
          nuclide_action.getParam<MooseFunctorName>("u");
      if (_problem->mesh().dimension() >= 2u)
        _moose_object_pars.set<MooseFunctorName>("v") =
            nuclide_action.getParam<MooseFunctorName>("v");
      if (_problem->mesh().dimension() >= 3u)
        _moose_object_pars.set<MooseFunctorName>("w") =
            nuclide_action.getParam<MooseFunctorName>("w");
    }

    _problem->addFVBC(_type, _name + "_" + nuclide + "_mass_fraction", _moose_object_pars);
  }
}
