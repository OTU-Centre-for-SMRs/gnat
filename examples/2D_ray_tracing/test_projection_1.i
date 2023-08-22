[Mesh]
  [NeutronicsDomain]
    type = FileMeshGenerator
    file = Storage_Radiation_Mesh_Tris.e
  []
[]

[AuxVariables]
  [UncollidedFlux]
    type = MooseVariable
    order = FIRST
    family = L2_LAGRANGE
  []
[]

[Transfers]
  [UncollidedFluxTransfer]
    type = MultiAppShapeEvaluationTransfer
    from_multi_app = Neutronics
    variable = UncollidedFlux
    source_variable = UncollidedFlux
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
