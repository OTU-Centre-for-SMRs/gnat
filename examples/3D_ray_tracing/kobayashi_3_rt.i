[Mesh]
  [domain]
    type = CartesianMeshGenerator
    dim = 3
    dx = '20 10 20 20 20 10 20'
    dy = '40 10 40 20 40 10 40'
    dz = '20 10 20 20 20 10 20'
    ix = '2 1 2 2 2 1 2'
    iy = '4 1 4 2 4 1 4'
    iz = '2 1 2 2 2 1 2'
    subdomain_id = '
    1 1 1 1 1 1 1
    1 1 1 1 1 1 1
    1 1 1 1 1 1 1
    1 1 1 1 1 1 1
    1 1 1 1 1 1 1
    1 1 1 1 1 1 1
    1 1 1 1 1 1 1

    1 0 1 1 1 0 1
    1 0 1 1 1 0 1
    1 1 1 1 1 1 1
    1 1 1 1 1 1 1
    1 1 1 1 1 1 1
    1 0 1 1 1 0 1
    1 0 1 1 1 0 1

    1 1 1 1 1 1 1
    1 0 1 1 1 0 1
    1 1 1 1 1 1 1
    1 1 1 1 1 1 1
    1 1 1 1 1 1 1
    1 0 1 1 1 0 1
    1 1 1 1 1 1 1

    1 1 1 1 1 1 1
    1 0 0 0 0 0 1
    1 1 1 0 1 1 1
    1 1 1 2 1 1 1
    1 1 1 0 1 1 1
    1 0 0 0 0 0 1
    1 1 1 1 1 1 1

    1 1 1 1 1 1 1
    1 0 1 1 1 0 1
    1 1 1 1 1 1 1
    1 1 1 1 1 1 1
    1 1 1 1 1 1 1
    1 0 1 1 1 0 1
    1 1 1 1 1 1 1

    1 0 1 1 1 0 1
    1 0 1 1 1 0 1
    1 1 1 1 1 1 1
    1 1 1 1 1 1 1
    1 1 1 1 1 1 1
    1 0 1 1 1 0 1
    1 0 1 1 1 0 1

    1 1 1 1 1 1 1
    1 1 1 1 1 1 1
    1 1 1 1 1 1 1
    1 1 1 1 1 1 1
    1 1 1 1 1 1 1
    1 1 1 1 1 1 1
    1 1 1 1 1 1 1'
  []
  uniform_refine = 1
[]

[UncollidedFlux]
  [Neutron]
    uncollided_flux_treatment = ray-tracing
    num_groups = 1
    max_anisotropy = 0

    volumetric_source_blocks = '2'
    volumetric_source_moments = '1.0'
    volumetric_source_anisotropies = '0'

    rt_n_polar = 30
    rt_n_azimuthal = 30
  []
[]

[TransportMaterials]
  #[Shield]
  #  type = ConstantNeutronicsMaterial
  #  transport_system = 'Neutron'
  #  group_absorption = '0.05'
  #  group_scattering = '0.05'
  #  anisotropy = 0
  #  group_speeds = '220000'
  #  block = '1 2'
  #[]
  #[Duct]
  #  type = ConstantNeutronicsMaterial
  #  transport_system = 'Neutron'
  #  group_absorption = '0.5e-4'
  #  group_scattering = '0.5e-4'
  #  anisotropy = 0
  #  group_speeds = '220000'
  #  block = '0'
  #[]
  [Shield]
    type = AbsorbingNeutronicsMaterial
    transport_system = 'Neutron'
    group_absorption = '0.1'
    group_speeds = '220000'
    block = '1 2'
  []
  [Duct]
    type = AbsorbingNeutronicsMaterial
    transport_system = 'Neutron'
    group_absorption = '1e-4'
    group_speeds = '220000'
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
