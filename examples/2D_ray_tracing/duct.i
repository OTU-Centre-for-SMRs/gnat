[Mesh]
  [domain]
    type = CartesianMeshGenerator
    dim = 2
    dx = '11 3 3 11'
    dy = '15 3 3 15'
    ix = '11 3 3 11'
    iy = '15 3 3 15'
    subdomain_id = '
    1 0 0 1
    1 2 2 1
    1 2 2 1
    1 0 0 1'
  []
  uniform_refine = 1
[]

[UncollidedFlux]
  [Neutron]
    uncollided_flux_treatment = ray-tracing
    num_groups = 2

    point_source_locations = '
    14.0  0.0 0.0
    14.0 36.0 0.0'
    point_source_moments = '
    1.0 0.0;
    1.0 0.0'
    point_source_anisotropies = '0 0'

    source_boundaries = '1 3'
    boundary_source_anisotropy = '0 0'
    boundary_source_moments = '0.5 0.1; 0.1 0.5'

    volumetric_source_blocks = '2'
    volumetric_source_anisotropies = '0'
    volumetric_source_moments = '1.0 0.5'

    rt_n_polar = 30
  []
[]

[TransportMaterials]
  [Shield]
    type = AbsorbingNeutronicsMaterial
    transport_system = 'Neutron'
    group_absorption = '0.5 1.0'
    group_speeds = '220000 220000'
    block = '1 2'
  []
  [Duct]
    type = AbsorbingNeutronicsMaterial
    transport_system = 'Neutron'
    group_absorption = '0.0 0.0'
    group_speeds = '220000 220000'
    block = '0'
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
[]
