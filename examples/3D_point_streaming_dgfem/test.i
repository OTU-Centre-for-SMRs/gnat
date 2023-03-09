# A simple test case with a purely absorbing medium and a point source in the
# middle of the domain.

[Mesh]
  [domain]
    type = CartesianMeshGenerator
    dim = 3
    dx = 10
    dy = 10
    dz = 10
    ix = 21
    iy = 21
    iz = 21
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

    point_source_locations = '5.0 5.0 5.0'
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
