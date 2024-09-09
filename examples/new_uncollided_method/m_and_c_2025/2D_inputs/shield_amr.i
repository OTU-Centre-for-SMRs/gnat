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
    group_total = 0.0
    source_location = '${src_x} ${src_y} ${src_z}'
  []
[]

[Materials]
  [Vacuum]
    type = AbsorbingTransportMaterial
    transport_system = ''
    num_groups = 1
    group_total = 0.0
    group_speeds = '1.0'
    block = '0'
    saaf_eta = 1.0
  []
  [Shield]
    type = AbsorbingTransportMaterial
    transport_system = ''
    num_groups = 1
    group_total = 1.0
    #group_total = 10.0
    group_speeds = '1.0'
    block = '1'
    saaf_eta = 1.0
  []
[]

[Adaptivity]
  marker = error_frac
  steps = 30

  [Indicators]
    #[error]
    #  type = GradientJumpIndicator
    #  variable = UncollidedFlux
    #[]
    #[error]
    #  type = SASFStreamingJumpIndicator
    #  variable = UncollidedFlux
    #  source_location = '${src_x} ${src_y} ${src_z}'
    #[]
    #[error]
    #  type = SASFConservationJumpIndicator
    #  variable = UncollidedFlux
    #  group_index = 0
    #  transport_system = ''
    #  source_location = '${src_x} ${src_y} ${src_z}'
    #[]
    [error]
      type = SASFTransverseJumpIndicator
      variable = UncollidedFlux
      source_location = '${src_x} ${src_y} ${src_z}'
    []
  []
  [Markers]
    [error_frac]
      type = ErrorFractionMarker
      coarsen = 0.0 # No coarsening for now.
      indicator = error
      #refine = 0.1
      #refine = 0.2
      #refine = 0.3
      #refine = 0.4
      refine = 0.5
    []
  []
[]

[Postprocessors/Num_Elements]
  type = NumElems
[]

[Executioner]
  type = Transient
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
  csv = true
[]
