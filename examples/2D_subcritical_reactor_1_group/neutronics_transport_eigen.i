[Mesh]
  [NeutronicsDomain]
    type = FileMeshGenerator
    file = Radiation_Mesh.e
  []
[]

[TransportSystems]
  [Neutron]
    num_groups = 1
    max_anisotropy = 0
    scheme = saaf_cfem
    particle_type = neutron

    flux_moment_names = 'Flux_Moment'
    debug_disable_fission = false
    eigen = true
    constant_ic = '1.0'

    order = FIRST
    family = LAGRANGE

    n_azimuthal = 5
    n_polar = 5

    vacuum_boundaries = 'vacuum'
  []
[]

[TransportMaterials]
  [Air]
    type = FileTransportMaterial
    transport_system = Neutron
    file_name = '../../data/mgxs/sc_1g_xs_macro.xml'
    source_material_id = '5'
    block = '2 7 8 9 ' #
  []

  [Cladding]
    type = FileTransportMaterial
    transport_system = Neutron
    file_name = '../../data/mgxs/sc_1g_xs_macro.xml'
    source_material_id = '7' #
    block = 4
  []

  [Fuel]
    type = FileTransportMaterial
    transport_system = Neutron
    file_name = '../../data/mgxs/sc_1g_xs_macro.xml'
    source_material_id = '6' #
    block = 3
  []

  [Graphite]
    type = FileTransportMaterial
    transport_system = Neutron
    file_name = '../../data/mgxs/sc_1g_xs_macro.xml'
    source_material_id = '8' #
    block = '10 11'
  []

  [Control]
    type = FileTransportMaterial
    transport_system = Neutron
    file_name = '../../data/mgxs/sc_1g_xs_macro.xml'
    source_material_id = '9' #
    block = 6
  []
[]

[Executioner]
  type = Eigenvalue

  initial_eigenvalue = 1.2

  solve_type = PJFNK
  petsc_options_iname = '-ksp_gmres_restart'
  petsc_options_value = ' 600'

  nl_abs_tol = 1e-12

  automatic_scaling = true
  off_diagonals_in_auto_scaling = true
  compute_scaling_once = false
[]

[Outputs]
  exodus = true
  execute_on = 'TIMESTEP_END'
[]
