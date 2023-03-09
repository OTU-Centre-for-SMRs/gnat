# A simple test case with a purely absorbing medium and a point source in the
# middle of the domain.

[Mesh]
  [domain]
    type = CartesianMeshGenerator
    dim = 2
    dx = 10
    dy = 10
    ix = 101
    iy = 101
  []
[]

[TransportSystems]
  [Neutron]
    num_groups = 1
    scheme = upwinding_dfem
    particle_type = neutron

    order = FIRST
    family = MONOMIAL

    n_azimuthal = 3
    n_polar = 3

    max_anisotropy = 0
    vacuum_boundaries = 'left right top bottom'

    point_source_locations = '5.0 5.0 0.0'
    point_source_intensities = '1000.0'
    point_source_groups = '1'

    debug_disable_scattering = true
  []
[]

[TransportMaterials]
  [Domain]
    type = AbsorbingNeutronicsMaterial
    transport_system = Neutron
    group_absorption = 2.0
    group_speeds = 2200.0
  []
[]

[Problem]
  type = FEProblem
[]

[Executioner]
  type = Steady
  solve_type = NEWTON
  nl_abs_tol = 1e-12
[]
