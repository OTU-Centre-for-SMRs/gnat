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
    scheme = upwinding_dfem
    particle_type = neutron
    num_groups = 1
    output_angular_fluxes = true

    order = FIRST
    family = MONOMIAL

    n_azimuthal = 1
    n_polar = 1

    max_anisotropy = 0
    vacuum_boundaries = 'right'
    reflective_boundaries = 'left'

    point_source_locations = '5.0 0.0 0.0'
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
    group_absorption = 0.0
    group_speeds = 2200.0
    #saaf_eta = 0.0
  []
[]

[Problem]
  type = FEProblem
[]

[Executioner]
  type = Steady
  solve_type = NEWTON
[]
