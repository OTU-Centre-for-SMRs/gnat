#include "UncollidedFluxAction.h"

#include "AddOutputAction.h"
#include "ActionWarehouse.h"

#include "RayKernelBase.h"

// Required for all treatments.
registerMooseAction("GnatApp", UncollidedFluxAction, "add_aux_variable");
registerMooseAction("GnatApp", UncollidedFluxAction, "add_aux_kernel");
registerMooseAction("GnatApp", UncollidedFluxAction, "add_postprocessor");

// Ray traced uncollided flux treatment.
registerMooseAction("GnatApp", UncollidedFluxAction, "add_user_object");
registerMooseAction("GnatApp", UncollidedFluxAction, "add_ray_kernel");

// SASF uncollided flux treatment.
registerMooseAction("GnatApp", UncollidedFluxAction, "add_variable");
registerMooseAction("GnatApp", UncollidedFluxAction, "add_kernel");
registerMooseAction("GnatApp", UncollidedFluxAction, "add_bc");

// For validating the problem.
registerMooseAction("GnatApp", UncollidedFluxAction, "post_mesh_prepared");

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
                             "Self-Adjoint Scalar Flux advection-reaction equation.");
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
  params.addParam<bool>("is_conservative_transfer_src",
                        true,
                        "Whether this transport action is providing flux moments to a sub/parent "
                        "app using conservative transfers. Setting this option to 'true' adds "
                        "post-processors to ensure flux moments are conservative.");
  //----------------------------------------------------------------------------
  // Boundary sources.
  params.addParam<std::vector<BoundaryName>>("source_boundaries",
                                             std::vector<BoundaryName>(),
                                             "The boundaries to apply incoming "
                                             "flux boundary conditions.");
  params.addParam<std::vector<std::vector<Real>>>(
      "boundary_source_moments",
      std::vector<std::vector<Real>>(),
      "A double vector containing the external source moments for "
      "all boundaries. The exterior vector must correspond with the surface source boundary "
      "conditions provided in 'source_boundaries'.");
  params.addParam<std::vector<unsigned int>>(
      "boundary_source_anisotropy",
      std::vector<unsigned int>(),
      "The degree of anisotropy of the boundary source moments. The exterior vector must "
      "correspond with the surface source boundary "
      "conditions provided in 'source_boundaries'.");
  params.addParamNamesToGroup(
      "source_boundaries boundary_source_moments boundary_source_anisotropy", "Boundary Sources");
  //----------------------------------------------------------------------------
  // Volumetric sources.
  params.addParam<std::vector<SubdomainName>>("volumetric_source_blocks",
                                              std::vector<SubdomainName>(),
                                              "The list of blocks (ids or "
                                              "names) that host a volumetric source.");
  params.addParam<std::vector<std::vector<Real>>>(
      "volumetric_source_moments",
      std::vector<std::vector<Real>>(),
      "A double vector containing a list of external source moments for all volumetric particle "
      "sources. The external vector should correspond with the order of "
      "'volumetric_source_blocks'.");
  params.addParam<std::vector<unsigned int>>(
      "volumetric_source_anisotropies",
      std::vector<unsigned int>(),
      "The anisotropies of the volumetric sources. The vector should correspond with the order of "
      "'volumetric_source_blocks'");
  params.addParamNamesToGroup(
      "volumetric_source_blocks volumetric_source_moments volumetric_source_anisotropies",
      "Volumetric Source");
  //----------------------------------------------------------------------------
  // Point sources.
  params.addParam<std::vector<Point>>("point_source_locations",
                                      std::vector<Point>(),
                                      "The locations of all isotropic "
                                      "point sources in the problem "
                                      "space.");
  params.addParam<std::vector<std::vector<Real>>>(
      "point_source_moments",
      std::vector<std::vector<Real>>(),
      "A double vector containing a list of external source moments for "
      "all point particle "
      "sources. The external vector should correspond with the order of "
      "'point_source_locations'.");
  params.addParam<std::vector<unsigned int>>(
      "point_source_anisotropies",
      std::vector<unsigned int>(),
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
                                            10,
                                            "rt_n_polar > 0",
                                            "Number of Legendre polar "
                                            "quadrature points in a single "
                                            "octant of the unit sphere. "
                                            "Defaults to 30.");
  params.addRangeCheckedParam<unsigned int>("rt_n_azimuthal",
                                            10,
                                            "rt_n_azimuthal > 0",
                                            "Number of Chebyshev azimuthal "
                                            "quadrature points in a single "
                                            "octant of the unit sphere. "
                                            "Defaults to 30.");
  params.addParamNamesToGroup("uncollided_flux_treatment rt_volume_order rt_face_order "
                              "rt_volume_type rt_face_type rt_n_polar rt_n_azimuthal",
                              "Ray Traced Uncollided Flux Treatment");

  //----------------------------------------------------------------------------
  // Parameters specific to the SASF approach.
  params.addParam<MooseEnum>("family",
                             AddVariableAction::getNonlinearVariableFamilies(),
                             "Specifies the family of FE shape functions to "
                             "use for this variable.");
  params.addParam<MooseEnum>("order",
                             AddVariableAction::getNonlinearVariableOrders(),
                             "Specifies the order of the FE shape "
                             "function to use for this variable "
                             "(additional orders not listed are "
                             "allowed).");
  params.addParam<Real>("scaling",
                        1.0,
                        "Specifies a scaling factor to apply to "
                        "this variable.");

  params.addParam<BoundaryName>("sasf_near_source_boundary",
                                "",
                                "The near-source boundary for mitigating singularities in "
                                "the Self-Adjoint Scalar Flux approach.");
  params.addParam<std::vector<BoundaryName>>("sasf_vacuum_boundaries",
                                             std::vector<BoundaryName>(),
                                             "The implicit vacuum boundaries in "
                                             "the Self-Adjoint Scalar Flux approach.");
  params.addParam<std::vector<Real>>(
      "sasf_near_source_cross_sections",
      std::vector<Real>(),
      "The groupwise total cross-sections in the near-source region.");

  return params;
}

UncollidedFluxAction::UncollidedFluxAction(const InputParameters & parameters)
  : GnatBaseAction(parameters),
    _num_groups(getParam<unsigned int>("num_groups")),
    _max_eval_anisotropy(getParam<unsigned int>("max_anisotropy")),
    _num_group_moments(0u),
    _uncollided_var_base_name(getParam<std::string>("uncollided_flux_moment_names")),
    _uncollided_treatment(
        getParam<MooseEnum>("uncollided_flux_treatment").getEnum<UncollidedTreatment>()),
    _conservative_src(getParam<bool>("is_conservative_transfer_src")),
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
        getParam<std::vector<unsigned int>>("volumetric_source_anisotropies")),
    _sasf_near_source_boundary(getParam<BoundaryName>("sasf_near_source_boundary")),
    _sasf_vacuum_boundaries(getParam<std::vector<BoundaryName>>("sasf_vacuum_boundaries")),
    _sasf_near_source_cross_sections(getParam<std::vector<Real>>("sasf_near_source_cross_sections"))
{
  if (_uncollided_treatment == UncollidedTreatment::SASF)
  {
    if (_point_source_locations.size() != 1u)
      mooseError("Only a single point source is supported when using the Self-Adjoint Scalar "
                 "Flux approach.");

    if (_boundary_source_moments.size() > 0u)
      mooseWarning("The Self-Adjoint Scalar Flux approach only supports point sources. All "
                   "boundary sources will be ignored.");

    if (_volumetric_source_blocks.size() > 0u)
      mooseWarning("The Self-Adjoint Scalar Flux approach only supports point sources. All "
                   "volume sources will be ignored.");

    if (_sasf_near_source_boundary == "")
      mooseError("A near source region boundary must be specified if using the Self-Adjoint Scalar "
                 "Flux approach.");

    if (_sasf_vacuum_boundaries.size() == 0u)
      mooseError("Vacuum boundaries must be provided if using the Self-Adjoint Scalar "
                 "Flux approach.");

    if (_sasf_near_source_cross_sections.size() != _num_groups)
      mooseError("The number of provided total cross sections for the Self-Adjoint Scalar Flux "
                 "approach must match the number of energy groups.");
  }
}

void
UncollidedFluxAction::act()
{
  if (_current_task == "post_mesh_prepared" && _uncollided_treatment == UncollidedTreatment::RT)
  {
    switch (_mesh->dimension())
    {
      case 2u:
        _num_group_moments = (_max_eval_anisotropy + 1u) * (_max_eval_anisotropy + 2u) / 2u;
        mooseWarning("The ray traced method only supports flatland calculations in 2D.");
        break;

      case 3u:
        _num_group_moments = (_max_eval_anisotropy + 1u) * (_max_eval_anisotropy + 1u);
        break;

      default:
        mooseError("Invalid mesh dimensionality (" + Moose::stringify(_mesh->dimension()) + ").");
        break;
    }
  }
  else
  {
    switch (_mesh->dimension())
    {
      case 2u:
        mooseWarning("The SASF method only supports flatland calculations in 2D.");
        break;

      case 3u:
        _num_group_moments = (_max_eval_anisotropy + 1u) * (_max_eval_anisotropy + 1u);
        break;

      default:
        mooseError("Invalid mesh dimensionality (" + Moose::stringify(_mesh->dimension()) + ").");
        break;
    }
  }

  switch (_uncollided_treatment)
  {
    case UncollidedTreatment::RT:
      actUncollidedFluxRT();
      break;

    case UncollidedTreatment::SASF:
      actUncollidedFluxSASF();
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

  if (_current_task == "add_postprocessor" && _conservative_src)
  {
    debugOutput("    - Adding the uncollided flux post-processors...");
    addUncollidedRayPostProcessors();
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
    params.set<unsigned int>("max_anisotropy") = _max_eval_anisotropy;
    params.set<unsigned int>("num_group_moments") = _num_group_moments;
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
    params.set<unsigned int>("components") = _num_groups * _num_group_moments;

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
      if (_mesh->dimension() == 2u)
      {
        for (unsigned int l = 0; l <= _max_eval_anisotropy; ++l)
        {
          for (int m = 0; m <= static_cast<int>(l); ++m)
          {
            const std::string var_name = _uncollided_var_base_name + "_" +
                                         Moose::stringify(g + 1u) + "_" + Moose::stringify(l) +
                                         "_" + Moose::stringify(m);

            _problem->addAuxVariable("MooseVariableConstMonomial", var_name, params);
            debugOutput("      - Adding auxvariable MooseVariableConstMonomial " + var_name + ".");
          }
        }
      }
      else
      {
        for (unsigned int l = 0; l <= _max_eval_anisotropy; ++l)
        {
          for (int m = -1 * static_cast<int>(l); m <= static_cast<int>(l); ++m)
          {
            const std::string var_name = _uncollided_var_base_name + "_" +
                                         Moose::stringify(g + 1u) + "_" + Moose::stringify(l) +
                                         "_" + Moose::stringify(m);

            _problem->addAuxVariable("MooseVariableConstMonomial", var_name, params);
            debugOutput("      - Adding auxvariable MooseVariableConstMonomial " + var_name + ".");
          }
        }
      }
    }
  } // MooseVariableConstMonomial
}

void
UncollidedFluxAction::addUncollidedRayAuxKernels()
{
  // Add ArrayVariableComponent.
  unsigned int index = 0u;
  for (unsigned int g = 0u; g < _num_groups; ++g)
  {
    if (_mesh->dimension() == 2u)
    {
      for (unsigned int l = 0; l <= _max_eval_anisotropy; ++l)
      {
        for (int m = 0; m <= static_cast<int>(l); ++m)
        {
          const std::string var_name = _uncollided_var_base_name + "_" + Moose::stringify(g + 1u) +
                                       "_" + Moose::stringify(l) + "_" + Moose::stringify(m);

          auto params = _factory.getValidParams("ArrayVariableComponent");
          params.set<std::vector<VariableName>>("array_variable")
              .emplace_back("RTUncollidedStorage");
          params.set<AuxVariableName>("variable") = var_name;
          params.set<unsigned int>("component") = index;

          if (isParamValid("block"))
          {
            params.set<std::vector<SubdomainName>>("block") =
                getParam<std::vector<SubdomainName>>("block");
          }

          params.set<ExecFlagEnum>("execute_on") = EXEC_TIMESTEP_END;

          _problem->addAuxKernel(
              "ArrayVariableComponent", "ArrayVariableComponent_" + var_name, params);
          debugOutput("      - Adding auxkernel ArrayVariableComponent for " + var_name + ".");

          index++;
        }
      }
    }
    else
    {
      for (unsigned int l = 0; l <= _max_eval_anisotropy; ++l)
      {
        for (int m = -1 * static_cast<int>(l); m <= static_cast<int>(l); ++m)
        {
          const std::string var_name = _uncollided_var_base_name + "_" + Moose::stringify(g + 1u) +
                                       "_" + Moose::stringify(l) + "_" + Moose::stringify(m);

          auto params = _factory.getValidParams("ArrayVariableComponent");
          params.set<std::vector<VariableName>>("array_variable")
              .emplace_back("RTUncollidedStorage");
          params.set<AuxVariableName>("variable") = var_name;
          params.set<unsigned int>("component") = index;

          if (isParamValid("block"))
          {
            params.set<std::vector<SubdomainName>>("block") =
                getParam<std::vector<SubdomainName>>("block");
          }

          params.set<ExecFlagEnum>("execute_on") = EXEC_TIMESTEP_END;

          _problem->addAuxKernel(
              "ArrayVariableComponent", "ArrayVariableComponent_" + var_name, params);
          debugOutput("      - Adding auxkernel ArrayVariableComponent for " + var_name + ".");

          index++;
        }
      }
    }
  } // ArrayVariableComponent
}

void
UncollidedFluxAction::addUncollidedRayPostProcessors()
{
  for (unsigned int g = 0u; g < _num_groups; ++g)
  {
    if (_mesh->dimension() == 2u)
    {
      for (unsigned int l = 0; l <= _max_eval_anisotropy; ++l)
      {
        for (int m = 0; m <= static_cast<int>(l); ++m)
        {
          const std::string var_name = _uncollided_var_base_name + "_" + Moose::stringify(g + 1u) +
                                       "_" + Moose::stringify(l) + "_" + Moose::stringify(m);

          auto params = _factory.getValidParams("ElementIntegralVariablePostprocessor");
          params.set<std::vector<VariableName>>("variable").emplace_back(var_name);
          params.set<ExecFlagEnum>("execute_on") = EXEC_TIMESTEP_END;
          params.set<std::vector<OutputName>>("outputs").emplace_back("none");

          if (isParamValid("block"))
            params.set<std::vector<SubdomainName>>("block") =
                getParam<std::vector<SubdomainName>>("block");

          _problem->addPostprocessor("ElementIntegralVariablePostprocessor",
                                     "ElementIntegralVariablePostprocessor_" + var_name + "_src",
                                     params);
          debugOutput(
              "      - Adding Transfer ElementIntegralVariablePostprocessor for the variable " +
              var_name + ".");
        }
      }
    }
    else
    {
      for (unsigned int l = 0; l <= _max_eval_anisotropy; ++l)
      {
        for (int m = -1 * static_cast<int>(l); m <= static_cast<int>(l); ++m)
        {
          const std::string var_name = _uncollided_var_base_name + "_" + Moose::stringify(g + 1u) +
                                       "_" + Moose::stringify(l) + "_" + Moose::stringify(m);

          auto params = _factory.getValidParams("ElementIntegralVariablePostprocessor");
          params.set<std::vector<VariableName>>("variable").emplace_back(var_name);
          params.set<ExecFlagEnum>("execute_on") = EXEC_TIMESTEP_END;
          params.set<std::vector<OutputName>>("outputs").emplace_back("none");

          if (isParamValid("block"))
            params.set<std::vector<SubdomainName>>("block") =
                getParam<std::vector<SubdomainName>>("block");

          _problem->addPostprocessor("ElementIntegralVariablePostprocessor",
                                     "ElementIntegralVariablePostprocessor_" + var_name + "_src",
                                     params);
          debugOutput(
              "      - Adding Transfer ElementIntegralVariablePostprocessor for the variable " +
              var_name + ".");
        }
      }
    }
  }
}

void
UncollidedFluxAction::actUncollidedFluxSASF()
{
  if (_current_task == "add_variable")
  {
    debugOutput("    - Adding the SASF uncollided flux variables...");
    addUncollidedSASFVariables();
  }

  if (_current_task == "add_kernel")
  {
    debugOutput("    - Adding the SASF uncollided flux kernels...");
    addUncollidedSASFKernels();
  }

  if (_current_task == "add_bc")
  {
    debugOutput("    - Adding the SASF uncollided flux BCs...");
    addUncollidedSASFBCs();
  }

  if (_current_task == "add_aux_variable" && _max_eval_anisotropy > 0)
  {
    debugOutput("    - Adding the SASF uncollided flux auxvariables...");
    addUncollidedSASFAuxVars();
  }

  if (_current_task == "add_aux_kernel" && _max_eval_anisotropy > 0)
  {
    debugOutput("    - Adding the SASF uncollided flux auxkernels...");
    addUncollidedSASFAuxKernels();
  }

  if (_current_task == "add_postprocessor" && _conservative_src)
  {
    debugOutput("    - Adding the SASF uncollided flux post-processors...");
    addUncollidedSASFPostProcessors();
  }
}

void
UncollidedFluxAction::addUncollidedSASFVariables()
{
  for (unsigned int g = 0u; g < _num_groups; ++g)
    addVariable(_uncollided_var_base_name + "_" + Moose::stringify(g + 1u) + "_0_0");
}

void
UncollidedFluxAction::addUncollidedSASFKernels()
{
  for (unsigned int g = 0u; g < _num_groups; ++g)
  {
    const std::string var_name =
        _uncollided_var_base_name + "_" + Moose::stringify(g + 1u) + "_0_0";

    // Add SASFAdvection.
    {
      auto params = _factory.getValidParams("SASFAdvection");
      params.set<NonlinearVariableName>("variable") = var_name;
      // Set the name of the TransportAction so it can fetch the appropriate material properties.
      params.set<std::string>("transport_system") = name();
      // Group index is required to fetch the group particle velocity.
      params.set<unsigned int>("group_index") = g;
      params.set<Point>("source_location") = _point_source_locations[0u];

      if (isParamValid("block"))
      {
        params.set<std::vector<SubdomainName>>("block") =
            getParam<std::vector<SubdomainName>>("block");
      }

      _problem->addKernel("SASFAdvection", "SASFAdvection_" + var_name, params);
      debugOutput("      - Adding kernel SASFAdvection for the variable " + var_name + ".");
    } // SASFAdvection

    // Add SASFRemoval.
    {
      auto params = _factory.getValidParams("SASFRemoval");
      params.set<NonlinearVariableName>("variable") = var_name;
      // Set the name of the TransportAction so it can fetch the appropriate material properties.
      params.set<std::string>("transport_system") = name();
      // Group index is required to fetch the group particle velocity.
      params.set<unsigned int>("group_index") = g;

      if (isParamValid("block"))
      {
        params.set<std::vector<SubdomainName>>("block") =
            getParam<std::vector<SubdomainName>>("block");
      }

      _problem->addKernel("SASFRemoval", "SASFRemoval_" + var_name, params);
      debugOutput("      - Adding kernel SASFRemoval for the variable " + var_name + ".");
    } // SASFRemoval
  }
}

void
UncollidedFluxAction::addUncollidedSASFBCs()
{
  for (unsigned int g = 0u; g < _num_groups; ++g)
  {
    const std::string var_name =
        _uncollided_var_base_name + "_" + Moose::stringify(g + 1u) + "_0_0";

    // Add SASFAnalyticalFluxBC.
    {
      auto params = _factory.getValidParams("SASFAnalyticalFluxBC");
      params.set<NonlinearVariableName>("variable") = var_name;
      // Group index is required to fetch the group particle velocity.
      params.set<unsigned int>("group_index") = g;
      params.set<std::vector<Real>>("group_total") = _sasf_near_source_cross_sections;
      params.set<std::vector<Real>>("group_source") = _point_source_moments[0u];
      params.set<Point>("source_location") = _point_source_locations[0u];

      params.set<std::vector<BoundaryName>>("boundary").emplace_back(_sasf_near_source_boundary);

      _problem->addBoundaryCondition(
          "SASFAnalyticalFluxBC", "SASFAnalyticalFluxBC_" + var_name, params);
      debugOutput("      - Adding BC SASFAnalyticalFluxBC for the variable " + var_name + ".");
    } // SASFAnalyticalFluxBC

    // Add SASFVacuumBC.
    {
      auto params = _factory.getValidParams("SASFVacuumBC");
      params.set<NonlinearVariableName>("variable") = var_name;
      params.set<Point>("source_location") = _point_source_locations[0u];

      params.set<std::vector<BoundaryName>>("boundary") = _sasf_vacuum_boundaries;

      _problem->addBoundaryCondition("SASFVacuumBC", "SASFVacuumBC_" + var_name, params);
      debugOutput("      - Adding BC SASFVacuumBC for the variable " + var_name + ".");
    } // SASFVacuumBC
  }
}

void
UncollidedFluxAction::addUncollidedSASFAuxVars()
{
  // Add AuxVariables.
  {
    auto fe_type = AddVariableAction::feType(_pars);
    auto type = AddVariableAction::variableType(fe_type, false, false);
    auto params = _factory.getValidParams(type);
    params.set<MooseEnum>("order") = fe_type.order.get_order();
    params.set<MooseEnum>("family") = Moose::stringify(fe_type.family);

    if (isParamValid("block"))
    {
      params.set<std::vector<SubdomainName>>("block") =
          getParam<std::vector<SubdomainName>>("block");
    }

    for (unsigned int g = 0; g < _num_groups; ++g)
    {
      if (_mesh->dimension() == 2u)
      {
        for (unsigned int l = 1; l <= _max_eval_anisotropy; ++l)
        {
          for (int m = 0; m <= static_cast<int>(l); ++m)
          {
            const std::string var_name = _uncollided_var_base_name + "_" +
                                         Moose::stringify(g + 1u) + "_" + Moose::stringify(l) +
                                         "_" + Moose::stringify(m);

            _problem->addAuxVariable(type, var_name, params);
            debugOutput("      - Adding auxvariable " + var_name + ".");
          }
        }
      }
      else
      {
        for (unsigned int l = 1; l <= _max_eval_anisotropy; ++l)
        {
          for (int m = -1 * static_cast<int>(l); m <= static_cast<int>(l); ++m)
          {
            const std::string var_name = _uncollided_var_base_name + "_" +
                                         Moose::stringify(g + 1u) + "_" + Moose::stringify(l) +
                                         "_" + Moose::stringify(m);

            _problem->addAuxVariable(type, var_name, params);
            debugOutput("      - Adding auxvariable " + var_name + ".");
          }
        }
      }
    }
  } // AuxVariables
}

void
UncollidedFluxAction::addUncollidedSASFAuxKernels()
{
  // Add ParticleFluxMomentSASF.
  for (unsigned int g = 0; g < _num_groups; ++g)
  {
    const auto scalar_flux = _uncollided_var_base_name + "_" + Moose::stringify(g + 1u) + "_0_0";
    if (_mesh->dimension() == 2u)
    {
      for (unsigned int l = 1; l <= _max_eval_anisotropy; ++l)
      {
        for (int m = 0; m <= static_cast<int>(l); ++m)
        {
          const std::string var_name = _uncollided_var_base_name + "_" + Moose::stringify(g + 1u) +
                                       "_" + Moose::stringify(l) + "_" + Moose::stringify(m);

          auto params = _factory.getValidParams("ParticleFluxMomentSASF");
          params.set<AuxVariableName>("variable") = var_name;
          params.set<Point>("source_location") = _point_source_locations[0u];
          params.set<std::vector<VariableName>>("uncollided_scalar_flux").emplace_back(scalar_flux);
          params.set<unsigned int>("degree") = l;
          params.set<int>("order") = m;

          _problem->addAuxKernel(
              "ParticleFluxMomentSASF", "ParticleFluxMomentSASF_" + var_name, params);
          debugOutput("      - Adding auxkernel ParticleFluxMomentSASF for the variable " +
                      var_name + ".");
        }
      }
    }
    else
    {
      for (unsigned int l = 1; l <= _max_eval_anisotropy; ++l)
      {
        for (int m = -1 * static_cast<int>(l); m <= static_cast<int>(l); ++m)
        {
          const std::string var_name = _uncollided_var_base_name + "_" + Moose::stringify(g + 1u) +
                                       "_" + Moose::stringify(l) + "_" + Moose::stringify(m);

          auto params = _factory.getValidParams("ParticleFluxMomentSASF");
          params.set<AuxVariableName>("variable") = var_name;
          params.set<Point>("source_location") = _point_source_locations[0u];
          params.set<std::vector<VariableName>>("uncollided_scalar_flux").emplace_back(scalar_flux);
          params.set<unsigned int>("degree") = l;
          params.set<int>("order") = m;

          _problem->addAuxKernel(
              "ParticleFluxMomentSASF", "ParticleFluxMomentSASF_" + var_name, params);
          debugOutput("      - Adding auxkernel ParticleFluxMomentSASF for the variable " +
                      var_name + ".");
        }
      }
    }
  } // ParticleFluxMomentSASF
}

void
UncollidedFluxAction::addUncollidedSASFPostProcessors()
{
  for (unsigned int g = 0u; g < _num_groups; ++g)
  {
    if (_mesh->dimension() == 2u)
    {
      for (unsigned int l = 0; l <= _max_eval_anisotropy; ++l)
      {
        for (int m = 0; m <= static_cast<int>(l); ++m)
        {
          const std::string var_name = _uncollided_var_base_name + "_" + Moose::stringify(g + 1u) +
                                       "_" + Moose::stringify(l) + "_" + Moose::stringify(m);

          auto params = _factory.getValidParams("ElementIntegralVariablePostprocessor");
          params.set<std::vector<VariableName>>("variable").emplace_back(var_name);
          params.set<ExecFlagEnum>("execute_on") = EXEC_TIMESTEP_END;
          params.set<std::vector<OutputName>>("outputs").emplace_back("none");

          if (isParamValid("block"))
            params.set<std::vector<SubdomainName>>("block") =
                getParam<std::vector<SubdomainName>>("block");

          _problem->addPostprocessor("ElementIntegralVariablePostprocessor",
                                     "ElementIntegralVariablePostprocessor_" + var_name + "_src",
                                     params);
          debugOutput(
              "      - Adding Transfer ElementIntegralVariablePostprocessor for the variable " +
              var_name + ".");
        }
      }
    }
    else
    {
      for (unsigned int l = 0; l <= _max_eval_anisotropy; ++l)
      {
        for (int m = -1 * static_cast<int>(l); m <= static_cast<int>(l); ++m)
        {
          const std::string var_name = _uncollided_var_base_name + "_" + Moose::stringify(g + 1u) +
                                       "_" + Moose::stringify(l) + "_" + Moose::stringify(m);

          auto params = _factory.getValidParams("ElementIntegralVariablePostprocessor");
          params.set<std::vector<VariableName>>("variable").emplace_back(var_name);
          params.set<ExecFlagEnum>("execute_on") = EXEC_TIMESTEP_END;
          params.set<std::vector<OutputName>>("outputs").emplace_back("none");

          if (isParamValid("block"))
            params.set<std::vector<SubdomainName>>("block") =
                getParam<std::vector<SubdomainName>>("block");

          _problem->addPostprocessor("ElementIntegralVariablePostprocessor",
                                     "ElementIntegralVariablePostprocessor_" + var_name + "_src",
                                     params);
          debugOutput(
              "      - Adding Transfer ElementIntegralVariablePostprocessor for the variable " +
              var_name + ".");
        }
      }
    }
  }
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
