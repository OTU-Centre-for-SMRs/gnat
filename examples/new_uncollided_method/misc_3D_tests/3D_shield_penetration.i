# A simple test case with a purely absorbing medium and a point source in the
# middle of the domain.

source = 1.0
xs = 0.0

src_x = 0.0
src_y = 0.0
src_z = 0.0

[GlobalParams]
  source_location = '${src_x} ${src_y} ${src_z}'
[]

[Mesh]
  [domain]
    type = FileMeshGenerator
    file = '3D_shield_penetration.e'
  []
  uniform_refine = 1
[]

[Adaptivity]
  marker = errorfrac
  steps = 15

  [Indicators]
    [error]
      type = GradientJumpIndicator
      variable = UncollidedFlux
      outputs = none
    []
  []

  [Markers]
    [errorfrac]
      type = ErrorFractionMarker
      indicator = error
      refine = 0.5
      coarsen = 0
      outputs = none
    []
  []
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
    group_source = ${source}
    group_total = ${xs}
  []
[]

[Materials]
  [Vacuum]
    type = AbsorbingNeutronicsMaterial
    transport_system = ''
    num_groups = 1
    group_total = ${xs}
    block = empty
  []

  [Shield]
    type = AbsorbingNeutronicsMaterial
    transport_system = ''
    num_groups = 1
    group_total = 1.0
    block = shield
  []
[]

[Executioner]
  type = Steady
  solve_type = PJFNK
  petsc_options_iname = '-ksp_gmres_restart'
  petsc_options_value = ' 600'

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
