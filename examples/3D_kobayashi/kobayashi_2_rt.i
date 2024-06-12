[Mesh]
  [domain]
    type = CartesianMeshGenerator
    dim = 3
    dx = '10 50'
    dy = '10 90'
    dz = '10 50'
    ix = '1 5'
    iy = '1 9'
    iz = '1 5'
    subdomain_id = '
    2 1
    0 1

    1 1
    1 1'
  []
  uniform_refine = 2
[]

[UncollidedFlux]
  [Neutron]
    uncollided_flux_treatment = ray-tracing
    num_groups = 1
    max_anisotropy = 0

    volumetric_source_blocks = '2'
    volumetric_source_moments = '1.0'
    volumetric_source_anisotropies = '0'

    rt_n_polar = 5
    rt_n_azimuthal = 5
  []
[]

[TransportMaterials]
  [Shield]
    type = ConstantTransportMaterial
    transport_system = 'Neutron'
    group_total = '0.1'
    group_scattering = '0.05'
    anisotropy = 0
    block = '1 2'
  []
  [Duct]
    type = ConstantTransportMaterial
    transport_system = 'Neutron'
    group_total = '1e-4'
    group_scattering = '0.5e-4'
    anisotropy = 0
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
