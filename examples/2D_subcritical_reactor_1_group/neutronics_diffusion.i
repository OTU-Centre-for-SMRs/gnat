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
    scheme = diffusion_cfem
    particle_type = neutron

    flux_moment_names = 'Flux_Moment'
    debug_disable_fission = false

    order = FIRST
    family = LAGRANGE

    vacuum_boundaries = 'vacuum'

    point_source_locations = '96.52 96.52 0.0'
    point_source_moments = '1e7'
    point_source_anisotropies = '0'
    scale_sources = true
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
  type = Steady
  solve_type = PJFNK
  petsc_options_iname = '-pc_type -pc_hypre_type -ksp_gmres_restart'
  petsc_options_value = ' hypre    boomeramg      300'

  l_max_its = 50
  nl_abs_tol = 1e-12

  automatic_scaling = true
  off_diagonals_in_auto_scaling = true
  compute_scaling_once = false
[]

[Outputs]
  exodus = true
[]
