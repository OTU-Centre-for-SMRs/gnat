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
    num_groups = 2
    max_anisotropy = 0
    scheme = saaf_cfem
    particle_type = neutron
    output_angular_fluxes = false

    order = FIRST
    family = LAGRANGE

    n_azimuthal = 2
    n_polar = 2

    vacuum_boundaries = 'left right top bottom'

    point_source_locations = '5.0 5.0 0.0'
    point_source_moments = '1e3 0.0'
    point_source_anisotropies = '0'
    scale_sources = true
  []
[]

[TransportMaterials]
  [Domain]
    type = ConstantTransportMaterial
    transport_system = Neutron
    anisotropy = 0
    group_total = '1.1 1.0'
    group_scattering = '0.5 0.5 0.0 1.0'
  []
[]

[Executioner]
  type = Steady
  solve_type = PJFNK
  petsc_options_iname = '-pc_type -pc_hypre_type -ksp_gmres_restart'
  petsc_options_value = ' hypre    boomeramg      10'
  l_max_its = 50
  nl_rel_tol = 1e-12
[]

[Outputs]
  exodus = true
  execute_on = 'TIMESTEP_END'
[]
