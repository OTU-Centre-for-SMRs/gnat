# A simple test case with a purely absorbing medium and a point source in the
# middle of the domain.

[Mesh]
  [domain]
    type = CartesianMeshGenerator
    dim = 1
    #dx = '1 9'
    #ix = '101 900'
    #subdomain_id = '0 1'
    dx = '10'
    ix = '1000'
  []
[]

[TransportSystems]
  [Neutron_1]
    num_groups = 1
    scheme = saaf_cfem
    particle_type = neutron
    output_angular_fluxes = false
    flux_moment_names = 'neutron_1'
    angular_flux_names = 'neutron_1'

    order = FIRST
    family = LAGRANGE

    n_polar = 10

    vacuum_boundaries = 'left right'

    point_source_locations = '5.0 0.0 0.0'
    point_source_moments = '1.0'
    point_source_anisotropies = '0'
  []
  [Neutron_2]
    num_groups = 1
    scheme = saaf_cfem
    particle_type = neutron
    output_angular_fluxes = false
    flux_moment_names = 'neutron_2'
    angular_flux_names = 'neutron_2'

    order = FIRST
    family = LAGRANGE

    n_polar = 20

    vacuum_boundaries = 'left right'

    point_source_locations = '5.0 0.0 0.0'
    point_source_moments = '1.0'
    point_source_anisotropies = '0'
  []
  [Neutron_3]
    num_groups = 1
    scheme = saaf_cfem
    particle_type = neutron
    output_angular_fluxes = false
    flux_moment_names = 'neutron_3'
    angular_flux_names = 'neutron_3'

    order = FIRST
    family = LAGRANGE

    n_polar = 40

    vacuum_boundaries = 'left right'

    point_source_locations = '5.0 0.0 0.0'
    point_source_moments = '1.0'
    point_source_anisotropies = '0'
  []
[]

[AuxVariables]
  [AnalyticalFluxField]
    type = MooseVariable
    order = FIRST
    family = LAGRANGE
  []
  [StoredSpatialError_1]
    type = MooseVariable
    order = FIRST
    family = LAGRANGE
  []
  [StoredSpatialError_2]
    type = MooseVariable
    order = FIRST
    family = LAGRANGE
  []
  [StoredSpatialError_3]
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
  [SpatialError_1]
    type = ParsedAux
    variable = StoredSpatialError_1
    coupled_variables = 'neutron_1_1_0_0 AnalyticalFluxField'
    expression = '100.0 * abs(neutron_1_1_0_0 - AnalyticalFluxField) / AnalyticalFluxField'
  []
  [SpatialError_2]
    type = ParsedAux
    variable = StoredSpatialError_2
    coupled_variables = 'neutron_2_1_0_0 AnalyticalFluxField'
    expression = '100.0 * abs(neutron_2_1_0_0 - AnalyticalFluxField) / AnalyticalFluxField'
  []
  [SpatialError_3]
    type = ParsedAux
    variable = StoredSpatialError_3
    coupled_variables = 'neutron_3_1_0_0 AnalyticalFluxField'
    expression = '100.0 * abs(neutron_3_1_0_0 - AnalyticalFluxField) / AnalyticalFluxField'
  []
[]

[Functions]
  [AnalyticalFlux]
    type = Analytical1DSlabFlux
    src_location = 5.0
    cross_section = 0.1
  []
[]

[TransportMaterials]
  [Domain_1]
    type = AbsorbingTransportMaterial
    transport_system = Neutron_1
    group_total = 0.1
  []
  [Domain_2]
    type = AbsorbingTransportMaterial
    transport_system = Neutron_2
    group_total = 0.1
  []
  [Domain_3]
    type = AbsorbingTransportMaterial
    transport_system = Neutron_3
    group_total = 0.1
  []
[]

[VectorPostprocessors]
  [LineAnal]
    type = LineFunctionSampler
    functions = AnalyticalFlux
    start_point = '5 0 0'
    end_point = '10 0 0'
    sort_by = x
    num_points = 100
    execute_on = timestep_end
  []
  [LineNum_1]
    type = LineValueSampler
    variable = neutron_1_1_0_0
    start_point = '5 0 0'
    end_point = '10 0 0'
    sort_by = x
    num_points = 100
    execute_on = timestep_end
  []
  [LineError_1]
    type = LineValueSampler
    variable = StoredSpatialError_1
    start_point = '5 0 0'
    end_point = '10 0 0'
    sort_by = x
    num_points = 100
    execute_on = timestep_end
  []
  [LineNum_2]
    type = LineValueSampler
    variable = neutron_2_1_0_0
    start_point = '5 0 0'
    end_point = '10 0 0'
    sort_by = x
    num_points = 100
    execute_on = timestep_end
  []
  [LineError_2]
    type = LineValueSampler
    variable = StoredSpatialError_2
    start_point = '5 0 0'
    end_point = '10 0 0'
    sort_by = x
    num_points = 100
    execute_on = timestep_end
  []
  [LineNum_3]
    type = LineValueSampler
    variable = neutron_3_1_0_0
    start_point = '5 0 0'
    end_point = '10 0 0'
    sort_by = x
    num_points = 100
    execute_on = timestep_end
  []
  [LineError_3]
    type = LineValueSampler
    variable = StoredSpatialError_3
    start_point = '5 0 0'
    end_point = '10 0 0'
    sort_by = x
    num_points = 100
    execute_on = timestep_end
  []
[]

[Outputs]
  exodus = true
  csv = true
  execute_on = TIMESTEP_END
[]

[Executioner]
  type = Steady
  solve_type = PJFNK
[]
