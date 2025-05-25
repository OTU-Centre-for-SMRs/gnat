[Mesh]
  [NeutronicsDomain]
    type = FileMeshGenerator
    file = 2D_quarter_core.e
  []
  uniform_refine = 0
[]

[TransportSystems]
  [Neutron]
    num_groups = 4
    max_anisotropy = 0
    scheme = saaf_cfem
    particle_type = neutron
    flux_moment_names = 'Flux_Moment'
    debug_disable_fission = false
    eigen = true

    order = FIRST
    family = LAGRANGE

    constant_ic = '1.0'

    n_azimuthal = 3
    n_polar = 3

    vacuum_boundaries = 'vacuum'
    reflective_boundaries = 'reflective'
  []
[]

[TransportMaterials]
  [Air]
    type = FileTransportMaterial
    transport_system = Neutron
    file_name = '../../data/mgxs/sc_4g_xs_macro.xml'
    source_material_id = '5' #
    block = 'air_gap src_air'
  []

  [Cladding]
    type = FileTransportMaterial
    transport_system = Neutron
    file_name = '../../data/mgxs/sc_4g_xs_macro.xml'
    source_material_id = '7' #
    block = 'cladding'
  []

  [Fuel]
    type = FileTransportMaterial
    transport_system = Neutron
    file_name = '../../data/mgxs/sc_4g_xs_macro.xml'
    source_material_id = '6' #
    block = 'fuel'
  []

  [Graphite]
    type = FileTransportMaterial
    transport_system = Neutron
    file_name = '../../data/mgxs/sc_4g_xs_macro.xml'
    source_material_id = '8' #
    block = 'core_graphite reflector_graphite'
  []

  [Control]
    type = FileTransportMaterial
    transport_system = Neutron
    file_name = '../../data/mgxs/sc_4g_xs_macro.xml'
    source_material_id = '8' #
    block = 'control_or_graphite'
  []
[]

[Executioner]
  type = Eigenvalue

  initial_eigenvalue = 1.2
  free_power_iterations = 6

  solve_type = PJFNK
  petsc_options_iname = '-pc_type -pc_hypre_type -ksp_gmres_restart'
  petsc_options_value = ' hypre    boomeramg      300'

  nl_abs_tol = 1e-12

  automatic_scaling = true
  off_diagonals_in_auto_scaling = true
  compute_scaling_once = false
[]

[VectorPostprocessors]
  [EigenValues]
    type = Eigenvalues
    inverse_eigenvalue = true
  []
[]

[Outputs]
  csv = true
  exodus = true
  execute_on = 'TIMESTEP_END'
[]
