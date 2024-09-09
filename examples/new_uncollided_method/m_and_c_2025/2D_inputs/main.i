[Mesh]
  [domain]
    type = FileMeshGenerator
    file = 'ray_traced_reference_sigma_1_out.e'
    #file = 'ray_traced_reference_sigma_10_out.e'
    use_for_exodus_restart = true
  []
[]

[AuxVariables]
  [SASFFlux]
    type = MooseVariable
    order = CONSTANT
    family = MONOMIAL
  []
  [RTFlux]
    type = MooseVariable
    order = CONSTANT
    family = MONOMIAL
    initial_from_file_var = uncollided_flux_moment_1_0_0
  []
  [UncollidedError]
    type = MooseVariable
    order = CONSTANT
    family = MONOMIAL
  []
  [UncollidedError_Abs]
    type = MooseVariable
    order = CONSTANT
    family = MONOMIAL
  []
[]

[AuxKernels]
  [SpatialError]
    type = ParsedAux
    variable = UncollidedError
    coupled_variables = 'SASFFlux RTFlux'
    expression = 'abs(SASFFlux - RTFlux) / RTFlux'
    block = '0 1'
  []

  [SpatialError_Abs]
    type = ParsedAux
    variable = UncollidedError_Abs
    coupled_variables = 'SASFFlux RTFlux'
    expression = 'abs(SASFFlux - RTFlux)'
    block = '0 1'
  []
[]

[Executioner]
  type = Transient
  num_steps = 30
[]

[Problem]
  type = FEProblem
  solve = false
  kernel_coverage_check = false
[]

[Transfers]
  [FetchSASFFlux]
    type = MultiAppGeneralFieldShapeEvaluationTransfer
    variable = SASFFlux
    from_multi_app = SASF
    source_variable = UncollidedFlux
    to_blocks = '0 1'
  []
[]

[Postprocessors]
  [Max_Rel_Error]
    type = ElementExtremeValue
    variable = UncollidedError
    value_type = max
    block = '0 1'
  []
  [Min_Rel_Error]
    type = ElementExtremeValue
    variable = UncollidedError
    value_type = min
    block = '0 1'
  []
  [Max_Abs_Error]
    type = ElementExtremeValue
    variable = UncollidedError_Abs
    value_type = max
    block = '0 1'
  []
  [Min_Abs_Error]
    type = ElementExtremeValue
    variable = UncollidedError_Abs
    value_type = min
    block = '0 1'
  []
  [Rel_Error_L2_Norm]
    type = ElementL2Norm
    variable = UncollidedError
    block = '0 1'
  []
  [Abs_Error_L2_Norm]
    type = ElementL2Norm
    variable = UncollidedError_Abs
    block = '0 1'
  []
[]

[MultiApps]
  [SASF]
    type = TransientMultiApp
    app_type = GnatApp
    input_files = 'shield_amr.i'
    execute_on = timestep_begin
  []
[]

[Outputs]
  exodus = true
  csv = true
[]

