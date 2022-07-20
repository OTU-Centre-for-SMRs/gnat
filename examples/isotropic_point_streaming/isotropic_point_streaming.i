# isotropic_point_streaming.i: A simple test case with a purely absorbing
# medium and a point source at the origin.

[Mesh]
  type = GeneratedMesh
  dim = 2
  nx = 10
  ny = 10 
  xmin = 0
  xmax = 1
  ymin = 0
  ymax = 1
[]

[NeutronActivationStudy]
  execution_type = steady
  family = LAGRANGE
  order = FIRST

  num_groups = 1
  max_anisotropy = 0
  vacuum_boundaries = '0 1 2 3'

  point_source_locations = '0.0 0.0 0.0'
  point_source_intensities = '1000.0'
  point_source_groups = '1'
[]

[Materials]
  [Domain]
    type = AbsorbingNeutronicsMaterial
    num_groups = 1
    group_removal = 1.0
    group_speeds = 1.0
  []
[]

[Problem]
  type = FEProblem # This is the "normal" type of Finite Element Problem in MOOSE
  coord_type = XYZ
[]

[Executioner]
  type = Steady # Steady state problem
  solve_type = NEWTON # Perform a Newton solve

  # Set PETSc parameters to optimize solver efficiency
  petsc_options_iname = '-pc_type -pc_hypre_type' # PETSc option pairs with values below
  petsc_options_value = ' hypre    boomeramg'
[]

[Outputs]
  exodus = true # Output Exodus format
[]
