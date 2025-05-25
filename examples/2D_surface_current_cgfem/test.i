[Mesh]
  [gen]
    type = GeneratedMeshGenerator
    dim = 2
    xmin = 0
    xmax = 10
    ymin = 0
    ymax = 10
    nx = 51
    ny = 51
  []
[]

[TransportSystems]
  [Neutron]
    scheme = saaf_cfem
    particle_type = neutron
    num_groups = 2
    max_anisotropy = 1

    flux_moment_names = 'Flux_Moment'

    order = FIRST
    family = LAGRANGE

    n_azimuthal = 5
    n_polar = 5

    vacuum_boundaries = 'bottom'
    current_boundaries = 'top left right'

    boundary_current_anisotropy = '1 1 1'
    boundary_currents = '1e3 0.0; 1e3 0.0; 1e3 0.0'
    scale_sources = true
  []
[]

[TransportMaterials]
  [Water]
    type = FileTransportMaterial
    transport_system = Neutron
    file_name = '../../data/mgxs/water_3g_xs_macro.xml'
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
