# A simple test case with a purely absorbing medium and a point source in the
# middle of the domain.

[Mesh]
  [domain]
    type = CartesianMeshGenerator
    dim = 1
    dx = 10
    ix = 1001
  []
[]

[TransportSystems]
  [Neutron]
    num_groups = 1
    scheme = saaf_cfem
    particle_type = neutron
    output_angular_fluxes = false

    order = FIRST
    family = LAGRANGE

    n_polar = 100

    vacuum_boundaries = 'left right'

    point_source_locations = '5.0 0.0 0.0'
    point_source_moments = '1.0'
    point_source_anisotropies = '1'
    scale_sources = true
  []
[]

[TransportMaterials]
  [Domain]
    type = AbsorbingNeutronicsMaterial
    transport_system = Neutron
    group_absorption = 0.1
    group_speeds = 2200.0
  []
[]

[Problem]
  type = FEProblem
[]

[Outputs]
  exodus = true
[]

[Executioner]
  type = Steady
  solve_type = PJFNK
[]
