[Mesh]
  [domain]
    type = CartesianMeshGenerator
    dim = 2
    dx = '9.0 2.0 3.5 1.0 4.5'
    dy = '9.0 2.0 3.5 1.0 4.5'
    ix = '90 20 35 10 45'
    iy = '90 20 35 10 45'
    subdomain_id = '0 0 0 0 0
                    0 2 0 0 0
                    0 0 0 0 0
                    0 0 0 1 0
                    0 0 0 0 0'
  []
  uniform_refine = 2
[]

[UncollidedFlux]
  [Neutron]
    uncollided_flux_treatment = ray-tracing
    num_groups = 1
    max_anisotropy = 0

    point_source_locations = '10.0 10.0 0.0'
    point_source_moments = '1.0'
    point_source_anisotropies = '0'
  []
[]

[TransportMaterials]
  [Shield]
    type = AbsorbingNeutronicsMaterial
    transport_system = 'Neutron'
    group_total = '1.0'
    block = '1'
  []
  [Empty]
    type = AbsorbingNeutronicsMaterial
    transport_system = 'Neutron'
    group_total = '0.0'
    block = '0 2'
  []
[]

[AuxVariables]
  [SASFFlux]
    type = MooseVariable
    order = FIRST
    family = LAGRANGE
  []
  [UncollidedError]
    type = MooseVariable
    order = CONSTANT
    family = MONOMIAL
  []
[]

[AuxKernels]
  [SpatialError]
    type = ParsedAux
    variable = UncollidedError
    coupled_variables = 'SASFFlux uncollided_flux_moment_1_0_0'
    expression = 'abs(SASFFlux - uncollided_flux_moment_1_0_0) / uncollided_flux_moment_1_0_0'
    block = '0 1'
  []
[]

[Problem]
  type = FEProblem
  solve = false
  kernel_coverage_check = false
[]

[Executioner]
  type = Steady
[]

[Outputs]
  exodus = true
  execute_on = 'TIMESTEP_END'
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

[MultiApps]
  [SASF]
    type = FullSolveMultiApp
    app_type = GnatApp
    input_files = 'nsr_10_percent_4_quad_shield.i'
    execute_on = INITIAL
    positions = '10.0 10.0 0.0'
  []
[]
