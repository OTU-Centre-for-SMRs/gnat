# A simple test case with a purely absorbing medium and a point source in the
# middle of the domain.

src_x = 0.0
src_y = 0.0
src_z = 0.0

shield_xs = 1e0
void_xs = 0.0

[GlobalParams]
  source_location = '${src_x} ${src_y} ${src_z}'
[]

[Mesh]
  [domain]
    type = FileMeshGenerator
    file = 'spatial_convergence_2_quad.e'
  []
[]

[AuxVariables]
  [AnalyticalFluxField]
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
  [CompAnalytical]
    type = FunctionAux
    variable = AnalyticalFluxField
    function = AnalyticalFlux
  []
  [SpatialError]
    type = ParsedAux
    variable = UncollidedError
    coupled_variables = 'UncollidedFlux AnalyticalFluxField'
    expression = 'abs(UncollidedFlux - AnalyticalFluxField) / AnalyticalFluxField'
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
    group_source = 1.0
    group_total = 0.0
  []
[]

[Materials]
  [Vacuum]
    type = AbsorbingNeutronicsMaterial
    transport_system = ''
    num_groups = 1
    group_total = '${void_xs}'
    block = 'empty'
    saaf_eta = 1.0
  []
  [Shield]
    type = AbsorbingNeutronicsMaterial
    transport_system = ''
    num_groups = 1
    group_total = '${shield_xs}'
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

[Functions]
  [AnalyticalFlux]
    type = AxisymmetricPointScalarFlux
    num_dims = 2
    src_strength = 1.0
    src_location = '0 0 0'
    axisymmetric_radii = '2.0 4.0 6.0 8.0 20.0'
    axisymmetric_cross_sections = '${void_xs} ${shield_xs} ${void_xs} ${shield_xs} ${void_xs}'
  []
[]

[Postprocessors]
  [Error]
    type = ElementL2Error
    variable = UncollidedFlux
    function = AnalyticalFlux
  []
  [H]
    type = AverageElementMinSize
  []
[]

[Outputs]
  exodus = true
  csv = true
  execute_on = TIMESTEP_END
[]
