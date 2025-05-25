[Mesh]
  [NeutronicsDomain]
    type = FileMeshGenerator
    file = Radiation_Mesh.e
  []
[]

[TransportSystems]
  [Neutron]
    num_groups = 1
    max_anisotropy = 1
    scheme = saaf_cfem
    particle_type = neutron

    output_angular_fluxes = false
    flux_moment_names = 'Flux_Moment'

    order = FIRST
    family = LAGRANGE

    n_azimuthal = 4
    n_polar = 4

    vacuum_boundaries = 'vacuum'

    is_conservative_transfer_src = true

    uncollided_from_multi_app = 'Uncollided'
  []
[]

[TransportMaterials]
  [Air]
    type = FileTransportMaterial
    transport_system = 'Neutron'
    file_name = '../../data/mgxs/sc_st_1g_xs_macro.xml'
    source_material_id = '5'
    block = air
  []
  [Cladding]
    type = FileTransportMaterial
    transport_system = 'Neutron'
    file_name = '../../data/mgxs/sc_st_1g_xs_macro.xml'
    source_material_id = '6'
    block = cladding
  []
  [Fuel]
    type = FileTransportMaterial
    transport_system = 'Neutron'
    file_name = '../../data/mgxs/sc_st_1g_xs_macro.xml'
    source_material_id = '7'
    block = fuel
  []
  [Wood]
    type = FileTransportMaterial
    transport_system = 'Neutron'
    file_name = '../../data/mgxs/sc_st_1g_xs_macro.xml'
    source_material_id = '8'
    block = wood
  []
[]

[Executioner]
  type = Steady
  solve_type = PJFNK
  petsc_options_iname = '-pc_type -pc_hypre_type -ksp_gmres_restart'
  petsc_options_value = ' hypre    boomeramg      10'

  l_max_its = 50
  nl_abs_tol = 1e-12

  automatic_scaling = true
  off_diagonals_in_auto_scaling = true
  compute_scaling_once = false
[]

[MultiApps]
  [Uncollided]
    type = FullSolveMultiApp
    app_type = GnatApp
    input_files = 'neutronics_uncollided.i'
    execute_on = INITIAL
  []
[]

[Outputs]
  exodus = true
  execute_on = 'TIMESTEP_END'
[]
