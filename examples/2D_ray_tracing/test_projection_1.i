[Mesh]
  [domain]
    type = CartesianMeshGenerator
    dim = 2
    dx = '11 3 3 11'
    dy = '15 3 3 15'
    ix = '11 3 3 11'
    iy = '15 3 3 15'
    subdomain_id = '
    1 0 0 1
    1 2 2 1
    1 2 2 1
    1 0 0 1'
  []
  uniform_refine = 1
[]

[AuxVariables]
  [UncollidedFlux_1]
    type = MooseVariable
    order = FIRST
    family = L2_LAGRANGE
  []

  [UncollidedFlux_2]
    type = MooseVariable
    order = FIRST
    family = L2_LAGRANGE
  []
[]

[Transfers]
  [UncollidedFluxTransfer_1]
    type = MultiAppShapeEvaluationTransfer
    from_multi_app = Neutronics
    variable = UncollidedFlux_1
    source_variable = uncollided_flux_moment_1_0_0
  []

  [UncollidedFluxTransfer_2]
    type = MultiAppShapeEvaluationTransfer
    from_multi_app = Neutronics
    variable = UncollidedFlux_2
    source_variable = uncollided_flux_moment_2_0_0
  []
[]

[MultiApps]
  [Neutronics]
    type = FullSolveMultiApp
    app_type = GnatApp
    input_files = 'test_projection_2.i'
    execute_on = TIMESTEP_BEGIN
  []
[]

[Problem]
  type = FEProblem
  solve = false
  kernel_coverage_check = false
  material_coverage_check = false
[]

[Executioner]
  type = Steady
[]

[Outputs]
  exodus = true
[]
