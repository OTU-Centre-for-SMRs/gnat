#include "UncollidedFluxAction.h"

#include "AddOutputAction.h"
#include "ActionWarehouse.h"

#include "RayKernelBase.h"

// Required for all treatments.
registerMooseAction("GnatApp", UncollidedFluxAction, "add_aux_variable");
registerMooseAction("GnatApp", UncollidedFluxAction, "add_aux_kernel");

// Raytracing uncollided flux treatment.
registerMooseAction("GnatApp", UncollidedFluxAction, "add_user_object");
registerMooseAction("GnatApp", UncollidedFluxAction, "add_ray_kernel");

InputParameters
UncollidedFluxAction::validParams()
{
  auto params = GnatBaseAction::validParams();
  params.addClassDescription(
      "An action which sets up an uncollided flux treatment using either the semi-analytical "
      "method of ray tracing, or the SASF-modified Shands-Hanopy-Morel technique.");

  //----------------------------------------------------------------------------
  // Simulation parameters.
  params.addParam<MooseEnum>("uncollided_flux_treatment",
                             MooseEnum("ray-tracing sasf", "ray-tracing"),
                             "The type of uncollided flux treatment to apply to the solution. "
                             "Options are a semi-analytical treatment using ray tracing or the "
                             "SASF-modified Shands-Hanopy-Morel technique.");
  params.addRequiredParam<unsigned int>("num_groups",
                                        "The number of spectral energy groups in the "
                                        "problem.");
  params.addParam<unsigned int>("max_anisotropy",
                                0u,
                                "The maximum degree of anisotropy to evaluate. "
                                "Defaults to 0 for isotropic scattering.");
  params.addParam<std::string>("uncollided_flux_moment_names",
                               "uncollided_flux_moment",
                               "Variable names for the moments of the angular "
                               "flux. The output format for the group flux "
                               "moments will be of the form "
                               "{uncollided_flux_moment_names}_g_l_m.");
  //----------------------------------------------------------------------------
  // Boundary sources.
  params.addParam<std::vector<BoundaryName>>("source_boundaries",
                                             "The boundaries to apply incoming "
                                             "flux boundary conditions.");
  params.addParam<std::vector<std::vector<Real>>>(
      "boundary_source_moments",
      "A double vector containing the external source moments for "
      "all boundaries. The exterior vector must correspond with the surface source boundary "
      "conditions provided in 'source_boundaries'.");
  params.addParam<std::vector<unsigned int>>(
      "boundary_source_anisotropy",
      "The degree of anisotropy of the boundary source moments. The exterior vector must "
      "correspond with the surface source boundary "
      "conditions provided in 'source_boundaries'.");
  params.addParamNamesToGroup(
      "source_boundaries boundary_source_moments boundary_source_anisotropy", "Boundary Sources");
  //----------------------------------------------------------------------------
  // Volumetric sources.
  params.addParam<std::vector<SubdomainName>>("volumetric_source_blocks",
                                              "The list of blocks (ids or "
                                              "names) that host a volumetric source.");
  params.addParam<std::vector<std::vector<Real>>>(
      "volumetric_source_moments",
      "A double vector containing a list of external source moments for all volumetric particle "
      "sources. The external vector should correspond with the order of "
      "'volumetric_source_blocks'.");
  params.addParam<std::vector<unsigned int>>(
      "volumetric_source_anisotropies",
      "The anisotropies of the volumetric sources. The vector should correspond with the order of "
      "'volumetric_source_blocks'");
  params.addParamNamesToGroup(
      "volumetric_source_blocks volumetric_source_moments volumetric_source_anisotropies",
      "Volumetric Source");
  //----------------------------------------------------------------------------
  // Point sources.
  params.addParam<std::vector<Point>>("point_source_locations",
                                      "The locations of all isotropic "
                                      "point sources in the problem "
                                      "space.");
  params.addParam<std::vector<std::vector<Real>>>(
      "point_source_moments",
      "A double vector containing a list of external source moments for all point particle "
      "sources. The external vector should correspond with the order of "
      "'point_source_locations'.");
  params.addParam<std::vector<unsigned int>>(
      "point_source_anisotropies",
      "The anisotropies of the point sources. The vector should correspond with the order of "
      "'point_source_locations'");
  params.addParamNamesToGroup("point_source_locations point_source_moments "
                              "point_source_anisotropies",
                              "Point Source");

  //----------------------------------------------------------------------------
  // Parameters specific to ray tracing.
  // Ray tracing spatial quadrature parameters.
  MooseEnum qorders("CONSTANT FIRST SECOND THIRD FOURTH FIFTH SIXTH SEVENTH EIGHTH NINTH TENTH "
                    "ELEVENTH TWELFTH THIRTEENTH FOURTEENTH FIFTEENTH SIXTEENTH SEVENTEENTH "
                    "EIGHTTEENTH NINTEENTH TWENTIETH",
                    "CONSTANT");
  params.addParam<MooseEnum>("rt_volume_order",
                             qorders,
                             "The volume quadrature rule order. For simplicity the same quadrature "
                             "order is used for volumetric sources and target elements.");
  params.addParam<MooseEnum>("rt_face_order", qorders, "The face quadrature rule order.");
  MooseEnum qtypes("GAUSS GRID", "GAUSS");
  params.addParam<MooseEnum>("rt_volume_type",
                             qtypes,
                             "The volume quadrature type. For simplicity the same quadrature type "
                             "is used for both volumetric sources and target elements.");
  params.addParam<MooseEnum>("rt_face_type", qtypes, "The face quadrature type.");
  // Ray tracing angular quadrature parameters.
  params.addRangeCheckedParam<unsigned int>("rt_n_polar",
                                            30,
                                            "rt_n_polar > 0",
                                            "Number of Legendre polar "
                                            "quadrature points in a single "
                                            "octant of the unit sphere. "
                                            "Defaults to 30.");
  params.addRangeCheckedParam<unsigned int>("rt_n_azimuthal",
                                            30,
                                            "rt_n_azimuthal > 0",
                                            "Number of Chebyshev azimuthal "
                                            "quadrature points in a single "
                                            "octant of the unit sphere. "
                                            "Defaults to 30.");
  params.addParamNamesToGroup("uncollided_flux_treatment rt_volume_order rt_face_order "
                              "rt_volume_type rt_face_type rt_n_polar rt_n_azimuthal",
                              "Uncollided Flux Treatment");

  return params;
}

UncollidedFluxAction::UncollidedFluxAction(const InputParameters & parameters)
  : GnatBaseAction(parameters),
    _num_groups(getParam<unsigned int>("num_groups")),
    _max_eval_anisotropy(getParam<unsigned int>("max_anisotropy")),
    _uncollided_var_base_name(getParam<std::string>("uncollided_flux_moment_names")),
    _uncollided_treatment(
        getParam<MooseEnum>("uncollided_flux_treatment").getEnum<UncollidedTreatment>()),
    _point_source_locations(getParam<std::vector<Point>>("point_source_locations")),
    _point_source_moments(getParam<std::vector<std::vector<Real>>>("point_source_moments")),
    _point_source_anisotropy(getParam<std::vector<unsigned int>>("point_source_anisotropies")),
    _source_side_sets(getParam<std::vector<BoundaryName>>("source_boundaries")),
    _boundary_source_moments(getParam<std::vector<std::vector<Real>>>("boundary_source_moments")),
    _boundary_source_anisotropy(getParam<std::vector<unsigned int>>("boundary_source_anisotropy")),
    _volumetric_source_blocks(getParam<std::vector<SubdomainName>>("volumetric_source_blocks")),
    _volumetric_source_moments(
        getParam<std::vector<std::vector<Real>>>("volumetric_source_moments")),
    _volumetric_source_anisotropy(
        getParam<std::vector<unsigned int>>("volumetric_source_anisotropies"))
{
  for (unsigned int g = 0; g < _num_groups; ++g)
  {
    // Set up variable names for the group flux moments.
    _group_flux_moments.emplace(g, std::vector<VariableName>());
    _group_flux_moments[g].emplace_back(_uncollided_var_base_name + "_" + Moose::stringify(g + 1u) +
                                        "_" + Moose::stringify(0u) + "_" + Moose::stringify(0u));
  }
}

void
UncollidedFluxAction::act()
{
  switch (_uncollided_treatment)
  {
    case UncollidedTreatment::RT:
      actUncollidedFluxRT();
      break;

    case UncollidedTreatment::SASF:
      mooseError("This approach has not been implemented yet.");
      break;

    default:
      break;
  }
}

void
UncollidedFluxAction::actUncollidedFluxRT()
{
  if (_current_task == "add_user_object")
  {
    debugOutput("    - Adding the uncollided flux ray study...");
    addUncollidedRayStudies();
  }

  if (_current_task == "add_ray_kernel")
  {
    debugOutput("    - Adding the uncollided flux ray kernels...");
    addUncollidedRayKernels();
  }

  if (_current_task == "add_aux_variable")
  {
    debugOutput("    - Adding the uncollided flux auxvariables...");
    addUncollidedRayAuxVars();

    debugOutput("    - Modifying output variables...");
    modifyRTOutputs();
  }

  if (_current_task == "add_aux_kernel")
  {
    debugOutput("    - Adding the uncollided flux auxkernels...");
    addUncollidedRayAuxKernels();
  }
}

void
UncollidedFluxAction::addUncollidedRayKernels()
{
  // Add UncollidedFluxRayKernel.
  {
    auto params = _factory.getValidParams("UncollidedFluxRayKernel");
    params.set<AuxVariableName>("variable") = "RTUncollidedStorage";
    params.set<unsigned int>("num_groups") = _num_groups;
    // We pretend that the uncollided flux actions are transport systems. For all intensive
    // purposes, they are.
    params.set<std::string>("transport_system") = name();

    // Query for the ray uncollided flux ray tracing study generated by this action.
    std::vector<UserObject *> uos;
    auto query = _problem->theWarehouse().query().condition<AttribSystem>("UserObject");
    query.condition<AttribName>("UncollidedFluxRayStudy_RTUncollidedStorage");
    query.queryInto(uos);
    params.set<RayTracingStudy *>("_ray_tracing_study") = dynamic_cast<RayTracingStudy *>(uos[0]);

    if (isParamValid("block"))
    {
      params.set<std::vector<SubdomainName>>("block") =
          getParam<std::vector<SubdomainName>>("block");
    }

    _problem->addObject<RayKernelBase>(
        "UncollidedFluxRayKernel", "UncollidedFluxRayKernel_RTUncollidedStorage", params);
    debugOutput(
        "      - Adding ray kernel UncollidedFluxRayKernel for the variable RTUncollidedStorage.");
  } // UncollidedFluxRayKernel
}
void
UncollidedFluxAction::addUncollidedRayStudies()
{
  // Add UncollidedFluxRayStudy.
  {
    auto params = _factory.getValidParams("UncollidedFluxRayStudy");
    params.set<unsigned int>("num_groups") = _num_groups;

    params.set<MooseEnum>("volume_order") = getParam<MooseEnum>("rt_volume_order");
    params.set<MooseEnum>("face_order") = getParam<MooseEnum>("rt_face_order");
    params.set<MooseEnum>("volume_type") = getParam<MooseEnum>("rt_volume_type");
    params.set<MooseEnum>("face_type") = getParam<MooseEnum>("rt_face_type");

    params.set<unsigned int>("n_polar") = getParam<unsigned int>("rt_n_polar");
    params.set<unsigned int>("n_azimuthal") = getParam<unsigned int>("rt_n_azimuthal");

    // Point sources.
    params.set<std::vector<Point>>("point_source_locations") = _point_source_locations;
    params.set<std::vector<std::vector<Real>>>("point_source_moments") = _point_source_moments;
    params.set<std::vector<unsigned int>>("point_source_anisotropies") = _point_source_anisotropy;

    // Surface sources.
    params.set<std::vector<BoundaryName>>("source_boundaries") = _source_side_sets;
    params.set<std::vector<std::vector<Real>>>("boundary_source_moments") =
        _boundary_source_moments;
    params.set<std::vector<unsigned int>>("boundary_source_anisotropy") =
        _boundary_source_anisotropy;

    // Volume sources.
    params.set<std::vector<SubdomainName>>("volumetric_source_blocks") = _volumetric_source_blocks;
    params.set<std::vector<std::vector<Real>>>("volumetric_source_moments") =
        _volumetric_source_moments;
    params.set<std::vector<unsigned int>>("volumetric_source_anisotropies") =
        _volumetric_source_anisotropy;

    _problem->addUserObject(
        "UncollidedFluxRayStudy", "UncollidedFluxRayStudy_RTUncollidedStorage", params);
    debugOutput(
        "      - Adding ray study UncollidedFluxRayStudy for the variable RTUncollidedStorage.");
  } // UncollidedFluxRayStudy
}
void
UncollidedFluxAction::addUncollidedRayAuxVars()
{
  // Add ArrayMooseVariable.
  {
    auto params = _factory.getValidParams("ArrayMooseVariable");
    params.set<MooseEnum>("order") = "CONSTANT";
    params.set<MooseEnum>("family") = "MONOMIAL";
    params.set<unsigned int>("components") = _num_groups;

    if (isParamValid("block"))
    {
      params.set<std::vector<SubdomainName>>("block") =
          getParam<std::vector<SubdomainName>>("block");
    }

    _problem->addAuxVariable("ArrayMooseVariable", "RTUncollidedStorage", params);
    debugOutput("      - Adding auxvariable ArrayMooseVariable RTUncollidedStorage.");
  } // ArrayMooseVariable

  // Add MooseVariableConstMonomial.
  {
    auto params = _factory.getValidParams("MooseVariableConstMonomial");

    if (isParamValid("block"))
    {
      params.set<std::vector<SubdomainName>>("block") =
          getParam<std::vector<SubdomainName>>("block");
    }

    for (unsigned int g = 0; g < _num_groups; ++g)
    {
      _problem->addAuxVariable("MooseVariableConstMonomial", _group_flux_moments[g][0u], params);
      debugOutput("      - Adding auxvariable MooseVariableConstMonomial " +
                  _group_flux_moments[g][0u] + ".");
    }
  } // MooseVariableConstMonomial
}

void
UncollidedFluxAction::addUncollidedRayAuxKernels()
{
  // Add ArrayVariableComponent.
  for (unsigned int g = 0u; g < _num_groups; ++g)
  {
    auto params = _factory.getValidParams("ArrayVariableComponent");
    params.set<std::vector<VariableName>>("array_variable").emplace_back("RTUncollidedStorage");
    params.set<AuxVariableName>("variable") = _group_flux_moments[g][0u];
    params.set<unsigned int>("component") = g;

    if (isParamValid("block"))
    {
      params.set<std::vector<SubdomainName>>("block") =
          getParam<std::vector<SubdomainName>>("block");
    }

    params.set<ExecFlagEnum>("execute_on") = EXEC_TIMESTEP_END;

    _problem->addAuxKernel(
        "ArrayVariableComponent", "ArrayVariableComponent_" + _group_flux_moments[g][0u], params);
    debugOutput("      - Adding auxkernel ArrayVariableComponent for " +
                _group_flux_moments[g][0u] + ".");
  } // ArrayVariableComponent
}

void
UncollidedFluxAction::modifyRTOutputs()
{
  // Fetch all AddOutputAction's from the action warehouse.
  const auto & output_actions = _app.actionWarehouse().getActionListByName("add_output");
  for (const auto & act : output_actions)
  {
    // Extract the Output action.
    AddOutputAction * action = dynamic_cast<AddOutputAction *>(act);
    if (!action)
      continue;

    InputParameters & output_params = action->getObjectParams();
    if (output_params.have_parameter<std::vector<VariableName>>("hide"))
      output_params.set<std::vector<VariableName>>("hide").emplace_back("RTUncollidedStorage");
  }
}
