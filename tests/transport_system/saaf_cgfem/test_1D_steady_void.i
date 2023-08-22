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
    vacuum_boundaries = 'left right'

    point_source_locations = '5.0 0.0 0.0'
    point_source_moments = '1000.0'
    point_source_anisotropies = '1'

    debug_verbosity = level0
    debug_disable_scattering = true
  []
[]

[TransportMaterials]
  [Domain]
    type = VoidNeutronicsMaterial
    transport_system = Neutron
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
  petsc_options_iname = '-pc_type -pc_hypre_type'
  petsc_options_value = 'hypre boomeramg'
[]
