# A simple test case with a purely absorbing medium and a point source in the
# middle of the domain.

[GlobalParams]
  source_location = '0.0 0.0 0.0'
[]

[Mesh]
  type = AnnularMesh
  nr = 5
  nt = 20
  rmin = 0
  rmax = 5
  growth_r = 1.0
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
    type = UncollidedPointSource
    variable = UncollidedFlux
    group_index = 0
    transport_system = ''
    group_source = 1.0
  []
[]

[Kernels]
  [UncollidedFluxAdvection]
    type = UncollidedPointAdvection
    variable = UncollidedFlux
    group_index = 0
    transport_system = ''
  []
  [UncollidedFluxRemoval]
    type = UncollidedRemoval
    variable = UncollidedFlux
    group_index = 0
    transport_system = ''
  []
[]

[BCs]
  [UncollidedFluxVacuumBC]
    type = UncollidedPointVacuumBC
    variable = UncollidedFlux
    boundary = '1'
  []
[]

[Materials]
  #[Shield]
  #  type = AbsorbingNeutronicsMaterial
  #  transport_system = ''
  #  num_groups = 1
  #  group_absorption = '1.0'
  #  group_speeds = '220000'
  #  block = '1'
  #[]
  [Vacuum]
    type = AbsorbingNeutronicsMaterial
    transport_system = ''
    num_groups = 1
    group_absorption = '0.0'
    group_speeds = '220000'
    #block = '1'
  []
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

  automatic_scaling = true
  off_diagonals_in_auto_scaling = true
  compute_scaling_once = false
[]
