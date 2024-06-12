[Mesh]
  [domain]
    type = CartesianMeshGenerator
    dim = 3
    dx = '0.9 0.2 0.9'
    dy = '0.9 0.2 0.9'
    dz = '0.9 0.2 0.9'
    ix = '9 2 9'
    iy = '9 2 9'
    iz = '9 2 9'
    subdomain_id = '
    0 0 0
    0 0 0
    0 0 0

    0 0 0
    0 1 0
    0 0 0

    0 0 0
    0 0 0
    0 0 0'
  []
  uniform_refine = 3
[]

[UncollidedFlux]
  [Neutron]
    uncollided_flux_treatment = ray-tracing
    num_groups = 1
    max_anisotropy = 0

    point_source_locations = '1.0 1.0 1.0'
    point_source_moments = '1.0'
    point_source_anisotropies = '0'
  []
[]

[TransportMaterials]
  [Shield]
    type = ConstantTransportMaterial
    transport_system = 'Neutron'
    group_total = '1e-2'
    group_scattering = '0.0'
    anisotropy = 0
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
[]
