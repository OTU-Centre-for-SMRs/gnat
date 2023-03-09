# A simple test case with a purely absorbing medium and a point source in the
# middle of the domain.

[Mesh]
  [domain]
    type = CartesianMeshGenerator
    dim = 2
    dx = '4 2 4'
    dy = '4 2 4'
    ix = '40 20 40'
    iy = '40 20 40'
    subdomain_id = '
      1 1 1
      1 2 1
      1 1 1'
  []
[]

[TransportSystems]
  [Neutron]
    scheme = upwinding_dfem
    particle_type = neutron
    num_groups = 1

    order = FIRST
    family = MONOMIAL

    n_azimuthal = 2
    n_polar = 2

    max_anisotropy = 0
    vacuum_boundaries = 'left right top bottom'
  []
[]

[TransportMaterials]
  [Domain]
    type = VoidNeutronicsMaterial
    transport_system = Neutron
    group_speeds = 220000.0
    block = '1'
  []
  [Source]
    type = SourceNeutronicsMaterial
    transport_system = Neutron
    group_absorption = 0.0
    source_anisotropy = 0
    group_source = 1000.0
    anisotropy = 0
    group_scattering = 0.0
    group_speeds = 220000.0
    block = '2'
  []
[]

[Problem]
  type = FEProblem
[]

[Executioner]
  type = Steady
  solve_type = NEWTON
[]
