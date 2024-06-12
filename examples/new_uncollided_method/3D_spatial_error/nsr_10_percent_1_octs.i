# A simple test case with a purely absorbing medium and a point source in the
# middle of the domain.

source = 1.0
xs = 0.0

[GlobalParams]
  source_location = '0.0 0.0 0.0'
[]

[Mesh]
  [domain]
    type = FileMeshGenerator
    file = 'nsr_10_one_octs_quads.e'
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

[AuxVariables]
  [UncollidedFluxAnal]
    type = MooseVariable
    order = FIRST
    family = LAGRANGE
  []
  [UncollidedError]
    type = MooseVariable
    order = FIRST
    family = LAGRANGE
  []
[]

[AuxKernels]
  [CompAnalyticalFlux]
    type = FunctionAux
    variable = UncollidedFluxAnal
    function = Exact
  []
  [SpatialError]
    type = ParsedAux
    variable = UncollidedError
    coupled_variables = 'UncollidedFlux UncollidedFluxAnal'
    expression = '100.0 * abs((UncollidedFlux - UncollidedFluxAnal) / UncollidedFluxAnal)'
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
  []
[]

[Functions]
  [Exact]
    type = ParsedFunction
    expression = '${source} * exp(-1.0 * ${xs} * sqrt(x * x + y * y + z * z)) / (4.0 * pi * (x * x + y * y + z * z))'
  []
[]

[Postprocessors]
  [L2Error]
    type = ElementL2Error
    function = Exact
    variable = UncollidedFlux
  []
  [h]
    type = AverageElementSize
  []
[]

[Executioner]
  type = Steady
  solve_type = PJFNK
  petsc_options_iname = '-pc_type -pc_hypre_type -ksp_gmres_restart'
  petsc_options_value = ' hypre    boomeramg      600'

  l_max_its = 50
  nl_abs_tol = 1e-12

  automatic_scaling = true
  off_diagonals_in_auto_scaling = true
  compute_scaling_once = false
[]

[Outputs]
  exodus = true
  csv = true
  execute_on = TIMESTEP_END
[]
