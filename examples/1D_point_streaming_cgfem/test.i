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
    scheme = saaf_cfem
    particle_type = neutron
    output_angular_fluxes = true

    order = FIRST
    family = LAGRANGE

    n_azimuthal = 1
    n_polar = 1

    vacuum_boundaries = 'left right'

    point_source_locations = '5.0 0.0 0.0'
    point_source_intensities = '1000.0'
    point_source_groups = '1'
  []
[]

[TransportMaterials]
  [Domain]
    type = AbsorbingNeutronicsMaterial
    transport_system = Neutron
    group_absorption = 0.0
    group_speeds = 2200.0
  []
[]

[Problem]
  type = FEProblem
[]

[Executioner]
  type = Steady
  solve_type = NEWTON
[]
