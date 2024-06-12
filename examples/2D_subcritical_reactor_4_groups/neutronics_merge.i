# An input deck which merges an analytically computed near-source
# region flux with the numerically computed SASF approach.

[Mesh]
  [NeutronicsDomain]
    type = FileMeshGenerator
    file = 2D_quarter_core.e
  []
  uniform_refine = 2
[]

[AuxVariables]
  [Flux_Moment_1_0_0_uncollided_NSR]
    type = MooseVariable
    order = FIRST
    family = L2_LAGRANGE
  []
  [Flux_Moment_1_0_0_uncollided_FSR]
    type = MooseVariable
    order = FIRST
    family = L2_LAGRANGE
  []

  [uncollided_flux_moment_1_0_0]
    type = MooseVariable
    order = FIRST
    family = L2_LAGRANGE
    initial_condition = '0.0'
  []
  [uncollided_flux_moment_2_0_0]
    type = MooseVariable
    order = FIRST
    family = L2_LAGRANGE
    initial_condition = '0.0'
  []
  [uncollided_flux_moment_3_0_0]
    type = MooseVariable
    order = FIRST
    family = L2_LAGRANGE
    initial_condition = '0.0'
  []
  [uncollided_flux_moment_4_0_0]
    type = MooseVariable
    order = FIRST
    family = L2_LAGRANGE
    initial_condition = '0.0'
  []
[]

[AuxKernels]
  [NSR_Analytical_Flux]
    type = FunctionAux
    variable = Flux_Moment_1_0_0_uncollided_NSR
    function = AnalyticalFlux
    block = src_air
    execute_on = INITIAL
  []

  [FSR_Sum]
    type = SumAux
    variable = uncollided_flux_moment_1_0_0
    values = 'Flux_Moment_1_0_0_uncollided_NSR Flux_Moment_1_0_0_uncollided_FSR'
    execute_on = TIMESTEP_END
  []
[]

[Functions]
  [AnalyticalFlux]
    type = AxisymmetricPointScalarFlux
    num_dims = 2
    src_strength = 1.5
    src_location = '0 0 0'
    axisymmetric_radii = '15.0'
    axisymmetric_cross_sections = '9.683916522303286e-05'
  []
[]

[Transfers]
  [Group_1]
    type = MultiAppGeneralFieldShapeEvaluationTransfer
    variable = Flux_Moment_1_0_0_uncollided_FSR
    source_variable = uncollided_flux_moment_1_0_0
    from_multi_app = Uncollided
  []
[]

[Problem]
  type = FEProblem
  solve = false
  kernel_coverage_check = false
[]

[Executioner]
  type = Steady
[]

[MultiApps]
  [Uncollided]
    type = FullSolveMultiApp
    app_type = GnatApp
    input_files = 'neutronics_uncollided.i'
    execute_on = INITIAL
  []
[]
