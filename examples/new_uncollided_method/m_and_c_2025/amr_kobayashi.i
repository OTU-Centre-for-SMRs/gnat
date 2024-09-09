src_x = 0.0
src_y = 0.0
src_z = 0.0

[Mesh]
  [domain]
    type = FileMeshGenerator
    file = 'mesh_in.e'
  []
  uniform_refine = 0
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
    source_location = '${src_x} ${src_y} ${src_z}'
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
    source_location = '${src_x} ${src_y} ${src_z}'
  []

  [UncollidedFluxSourceBC]
    type = SASFAnalyticalFluxBC
    variable = UncollidedFlux
    boundary = 'nsr'
    group_index = 0
    group_source = 1.0
    group_total = 1e-4
    source_location = '${src_x} ${src_y} ${src_z}'
  []
[]

[Materials]
  [Duct]
    type = AbsorbingTransportMaterial
    transport_system = ''
    num_groups = 1
    group_total = 1e-4
    group_speeds = '1.0'
    block = '1'
    saaf_eta = 1.0
  []
  [Shield]
    type = AbsorbingTransportMaterial
    transport_system = ''
    num_groups = 1
    group_total = 1e-1
    group_speeds = '1.0'
    block = '0'
    saaf_eta = 1.0
  []
[]

[Adaptivity]
  marker = combo
  steps = 30
  max_h_level = 7

  [Indicators]
    #[error]
    #  type = SASFConservationJumpIndicator
    #  variable = UncollidedFlux
    #  group_index = 0
    #  transport_system = ''
    #  source_location = '${src_x} ${src_y} ${src_z}'
    #[]
    #[error]
    #  type = SASFStreamingJumpIndicator
    #  variable = UncollidedFlux
    #  source_location = '${src_x} ${src_y} ${src_z}'
    #[]
    [error]
      type = GradientJumpIndicator
      variable = UncollidedFlux
    []
  []
  [Markers]
    [error_frac]
      type = ErrorFractionMarker
      indicator = error
      refine = 0.3
      coarsen = 0.0 # No coarsening for now.
    []
    [negative]
      type = ValueThresholdMarker
      variable = UncollidedFlux
      refine = 0.0
      invert = true
    []
    [combo]
      type = ComboMarker
      markers = 'error_frac negative'
    []
    #[error_frac]
    #  type = ErrorFractionMarker
    #  indicator = error
    #  refine = 0.3
    #  coarsen = 0.0 # No coarsening for now.
    #[]
  []
[]

[Postprocessors/Num_Elements]
  type = NumElems
[]

[Executioner]
  type = Steady
  solve_type = PJFNK
  petsc_options_iname = '-ksp_gmres_restart'
  petsc_options_value = ' 600'

  l_max_its = 50
  #nl_abs_tol = 1e-12
  nl_rel_tol = 1e-8

  automatic_scaling = true
  off_diagonals_in_auto_scaling = true
  compute_scaling_once = false
[]

[Outputs]
  exodus = true
  csv = true
[]
