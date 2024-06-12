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
    file = 'nsr_10_percent_2_quads.e'
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

#[DiracKernels]
#  [UncollidedFluxPointSource]
#    type = SASFPointSource
#    variable = UncollidedFlux
#    group_index = 0
#    transport_system = ''
#    group_source = 1.0
#  []
#[]

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
    expression = '${source} * exp(-1.0 * ${xs} * sqrt(x * x + y * y)) / (2.0 * pi * sqrt(x * x + y * y))'
  []
[]

[Postprocessors]
  [error]
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
  execute_on = TIMESTEP_END
[]
