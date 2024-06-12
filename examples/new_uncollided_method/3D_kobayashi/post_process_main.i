[Mesh]
  [Domain]
    type = FileMeshGenerator
    file = 'kobayashi_3_rt_out.e'
    use_for_exodus_restart = true
  []
[]

[AuxVariables]
  [UncollidedRayTraced]
    type = MooseVariable
    order = CONSTANT
    family = MONOMIAL
    initial_from_file_var = 'uncollided_flux_moment_1_0_0'
    initial_from_file_timestep = LATEST
    block = '0 1'
  []
  [UncollidedSASF]
    type = MooseVariable
    order = FIRST
    family = LAGRANGE
    block = '0 1'
  []
  [UncollidedError]
    type = MooseVariable
    order = CONSTANT
    family = MONOMIAL
    block = '0 1'
  []
[]

[AuxKernels]
  [SpatialError]
    type = ParsedAux
    variable = UncollidedError
    coupled_variables = 'UncollidedSASF UncollidedRayTraced'
    expression = 'abs(UncollidedSASF - UncollidedRayTraced) / UncollidedRayTraced'
    block = '0 1'
  []
[]

[Transfers]
  [FetchSASFFlux]
    type = MultiAppGeneralFieldShapeEvaluationTransfer
    variable = UncollidedSASF
    from_multi_app = Sub
    source_variable = UncollidedFlux
    to_blocks = '0 1'
  []
[]

[MultiApps]
  [Sub]
    type = FullSolveMultiApp
    app_type = GnatApp
    input_files = 'post_process_sub.i'
    execute_on = INITIAL
  []
[]

[Problem]
  solve = false
  kernel_coverage_check = false
  material_coverage_check = false
[]

[Executioner]
  type = Steady
[]

[Outputs]
  exodus = true
  execute_on = TIMESTEP_END
[]
