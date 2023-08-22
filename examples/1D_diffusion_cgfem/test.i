# A simple test case with a purely absorbing medium and a point source in the
# middle of the domain.

[Mesh]
  [domain]
    type = CartesianMeshGenerator
    dim = 1
    dx = 10
    ix = 101
  []
[]

[TransportSystems]
  [Neutron]
    num_groups = 1
    scheme = diffusion_cfem
    particle_type = neutron

    order = FIRST
    family = LAGRANGE

    vacuum_boundaries = 'left right'

    point_source_locations = '5.0 0.0 0.0'
    point_source_moments = '1000.0'
    point_source_anisotropies = '0'
    scale_sources = true
  []
[]

[TransportMaterials]
  [Domain]
    type = AbsorbingNeutronicsMaterial
    transport_system = Neutron
    group_absorption = 0.5
    group_speeds = 2200.0
  []
[]

[Outputs]
  exodus = true
[]

[Executioner]
  type = Steady
  solve_type = PJFNK
[]
