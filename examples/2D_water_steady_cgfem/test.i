# A simple test case with a purely absorbing medium and a point source in the
# middle of the domain.

[Mesh]
  [domain]
    type = CartesianMeshGenerator
    dim = 2
    dx = '3.5 3.1 3.5'
    dy = '3.5 3.1 3.5'
    ix = '35 31 35'
    iy = '35 31 35'
    subdomain_id = '
      1 1 1
      1 2 1
      1 1 1'
  []
[]

[TransportSystems]
  [Neutron]
    num_groups = 2
    max_anisotropy = 1
    scheme = saaf_cfem
    particle_type = neutron

    order = FIRST
    family = LAGRANGE

    n_azimuthal = 5
    n_polar = 5

    vacuum_boundaries = 'left right top bottom'

    point_source_locations = '5.0 5.0 0.0'
    point_source_moments = '1e6 0.0'
    point_source_anisotropies = '0'
    scale_sources = true
  []
[]

[TransportMaterials]
  [Water]
    type = FileTransportMaterial
    transport_system = Neutron
    file_name = 'macro_xs.xml'
    source_material_id = '1'
  []
[]

[Executioner]
  type = Steady
  solve_type = PJFNK
  petsc_options_iname = '-pc_type -pc_hypre_type -ksp_gmres_restart'
  petsc_options_value = ' hypre    boomeramg      100'
  l_max_its = 50
  nl_rel_tol = 1e-12

  automatic_scaling = true
  off_diagonals_in_auto_scaling = true
  compute_scaling_once = false
[]

[Outputs]
  exodus = true
  execute_on = 'TIMESTEP_END'
[]
