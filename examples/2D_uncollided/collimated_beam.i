# A simple test case with a purely absorbing medium and a point source in the
# middle of the domain.

[GlobalParams]
  source_location = '0.0 3.0 0.0'
[]

[Mesh]
  [domain]
    type = FileMeshGenerator
    file = collimated_beam.e
  []
  second_order = false
  parallel_type = DISTRIBUTED
[]

[Variables]
  [UncollidedFlux]
    type = MooseVariable
    order = FIRST
    family = LAGRANGE
    initial_condition = 0.0
  []
[]

[DiracKernels]
  [UncollidedFluxPointSource]
    type = SASFPointSource
    variable = UncollidedFlux
    group_index = 0
    transport_system = ''
    group_source = 1.0
  []
[]

[Kernels]
  [UncollidedFluxAdvection]
    type = SASFAdvection
    variable = UncollidedFlux
    group_index = 0
    transport_system = ''
  []
  [UncollidedFluxRemoval]
    type = SASFRemoval
    variable = UncollidedFlux
    group_index = 0
    transport_system = ''
  []
[]

[BCs]
  [UncollidedFluxVacuumBC]
    type = SASFVacuumBC
    variable = UncollidedFlux
    boundary = '0 1 2 3'
  []
[]

[Materials]
  [Duct]
    type = AbsorbingNeutronicsMaterial
    transport_system = ''
    num_groups = 1
    group_absorption = '1e-6'
    group_speeds = '220000'
    block = '0'
    saaf_eta = 2.0
    saaf_c = 1.0
  []
  [Shield]
    type = AbsorbingNeutronicsMaterial
    transport_system = ''
    num_groups = 1
    group_absorption = '2e0'
    group_speeds = '220000'
    block = '1'
    saaf_eta = 2.0
    saaf_c = 1.0
  []
  #[Background]
  #  type = AbsorbingNeutronicsMaterial
  #  transport_system = ''
  #  num_groups = 1
  #  group_absorption = '1e-1'
  #  group_speeds = '220000'
  #  block = '2'
  #  saaf_eta = 2.0
  #  saaf_c = 1.0
  #[]
  #[Other]
  #  type = AbsorbingNeutronicsMaterial
  #  transport_system = ''
  #  num_groups = 1
  #  group_absorption = '0.5'
  #  group_speeds = '220000'
  #  block = '0'
  #[]
[]

[Problem]
  type = FEProblem
[]

[Outputs]
  exodus = true
[]

[Executioner]
  type = Steady
  solve_type = PJFNK
  #petsc_options_iname = '-pc_type -pc_hypre_type -ksp_gmres_restart'
  #petsc_options_value = 'hypre boomeramg 10'
  petsc_options_iname = '-pc_type -pc_factor_shift_type'
  petsc_options_value = 'lu       NONZERO'

  l_max_its = 50
  nl_abs_tol = 1e-12

  #automatic_scaling = true
  #off_diagonals_in_auto_scaling = true
  #compute_scaling_once = false
[]
