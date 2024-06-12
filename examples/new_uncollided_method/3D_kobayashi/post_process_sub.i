[Mesh]
  [Domain]
    type = FileMeshGenerator
    file = 'kobayashi_3_out.e'
    use_for_exodus_restart = true
  []
[]

[AuxVariables]
  [UncollidedFlux]
    type = MooseVariable
    order = FIRST
    family = LAGRANGE
    initial_from_file_var = 'UncollidedFlux'
    initial_from_file_timestep = LATEST
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
  exodus = false
[]
