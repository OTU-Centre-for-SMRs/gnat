[Mesh]
  [domain]
    type = CartesianMeshGenerator
    dim = 2
    dx = 10
    dy = 10
    ix = 10
    iy = 10
  []
[]

[TransportSystems]
  [Neutron]
    scheme = upwinding_dfem
    particle_type = neutron
    num_groups = 1
    output_angular_fluxes = true

    order = FIRST
    family = MONOMIAL

    n_azimuthal = 1
    n_polar = 1

    max_anisotropy = 0
    vacuum_boundaries = 'left right top bottom'

    point_source_locations = '5.0 5.0 0.0'
    point_source_intensities = '1000.0'
    point_source_groups = '1'

    debug_verbosity = level0
    debug_disable_scattering = true
  []
[]

[TransportMaterials]
  [Domain]
    type = AbsorbingNeutronicsMaterial
    transport_system = Neutron
    group_absorption = 1.0
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
