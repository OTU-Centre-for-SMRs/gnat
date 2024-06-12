# A simple test case with a purely absorbing medium and a point source in the
# middle of the domain.

[Mesh]
  [domain]
    type = CartesianMeshGenerator
    dim = 1
    dx = 10
    ix = 100
  []
[]

[TransportSystems]
  [Neutron]
    scheme = saaf_cfem
    particle_type = neutron
    num_groups = 1
    output_angular_fluxes = true

    order = FIRST
    family = LAGRANGE

    n_azimuthal = 1
    n_polar = 1

    max_anisotropy = 0
    vacuum_boundaries = 'right'
    reflective_boundaries = 'left'

    point_source_locations = '5.0 0.0 0.0'
    point_source_moments = '1000.0'
    point_source_anisotropies = '0'

    debug_verbosity = level0
    debug_disable_scattering = true
  []
[]

[TransportMaterials]
  [Domain]
    type = AbsorbingTransportMaterial
    transport_system = Neutron
    group_total = 0.0
    group_speeds = 2200.0
    #saaf_eta = 0.0
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
  solve_type = NEWTON
[]
