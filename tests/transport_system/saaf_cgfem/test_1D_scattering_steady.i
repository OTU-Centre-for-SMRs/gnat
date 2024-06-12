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
    use_scattering_jacobians = true

    order = FIRST
    family = LAGRANGE

    n_azimuthal = 1
    n_polar = 1

    max_anisotropy = 0
    vacuum_boundaries = 'left right'

    point_source_locations = '5.0 0.0 0.0'
    point_source_moments = '1000.0'
    point_source_anisotropies = '0'

    debug_verbosity = level0
  []
[]

[TransportMaterials]
  [Domain1]
    type = ConstantTransportMaterial
    transport_system = Neutron
    anisotropy = 0
    group_total = 2.0
    group_scattering = 1.0
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
  line_search = default
  petsc_options_iname = '-pc_type -pc_hypre_type -ksp_gmres_restart'
  petsc_options_value = 'hypre boomeramg 10'
  l_max_its = 50
  nl_rel_tol = 1e-12
[]
