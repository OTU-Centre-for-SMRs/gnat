# A simple test case with a purely absorbing medium and a point source in the
# middle of the domain.

src_x = 0.0
src_y = 0.0
src_z = 0.0

[GlobalParams]
  source_location = '${src_x} ${src_y} ${src_z}'
[]

[Mesh]
  [domain]
    type = FileMeshGenerator
    file = 'shield.e'
  []
  uniform_refine = 3
[]

[Variables]
  [UncollidedFlux]
    type = MooseVariable
    order = FIRST
    family = LAGRANGE
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
    boundary = 'vacuum'
  []

  [UncollidedFluxSourceBC]
    type = SASFAnalyticalFluxBC
    variable = UncollidedFlux
    boundary = 'analytical_flux'
    group_index = 0
    group_source = 1.0
    group_total = 0.0
  []
[]

[Materials]
  [Vacuum]
    type = AbsorbingNeutronicsMaterial
    transport_system = ''
    num_groups = 1
    group_total = 0.0
    block = 'empty'
    saaf_eta = 1.0
  []
  [Shield]
    type = AbsorbingNeutronicsMaterial
    transport_system = ''
    num_groups = 1
    group_total = 1.0
    block = 'shield'
    saaf_eta = 1.0
  []
[]

[Executioner]
  type = Steady
  solve_type = PJFNK
  petsc_options_iname = '-pc_type -pc_factor_shift_type'
  petsc_options_value = ' lu       NONZERO'

  l_max_its = 50
  nl_abs_tol = 1e-12

  automatic_scaling = true
  off_diagonals_in_auto_scaling = true
  compute_scaling_once = false
[]

[Outputs]
  exodus = true
  execute_on = TIMESTEP_END
[]
