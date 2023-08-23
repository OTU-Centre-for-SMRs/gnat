[Mesh]
  [NeutronicsDomain]
    type = FileMeshGenerator
    file = Full_Neutronics_3_Simple.e
  []
  uniform_refine = 2
[]

[TransportSystems]
  [Neutron]
    num_groups = 8
    max_anisotropy = 0
    scheme = saaf_cfem
    particle_type = neutron

    flux_moment_names = 'Flux_Moment'
    scale_sources = true

    order = FIRST
    family = LAGRANGE

    n_azimuthal = 5
    n_polar = 5

    vacuum_boundaries = 'vacuum'
    current_boundaries = 'neutron'

    boundary_current_anisotropy = '1'
    boundary_currents = '
    7.534300e7 1.060700e7 3.250000e5 6.347148e9
    3.566577e9 0.0        0.0        0.0'
  []
[]

[TransportMaterials]
  [RPV]
    type = FileNeutronicsMaterial
    transport_system = Neutron
    file_name = 'data/xs_macro.xml'
    source_material_id = '6'
    block = rpv
  []
  [Air]
    type = FileNeutronicsMaterial
    transport_system = Neutron
    file_name = 'data/xs_macro.xml'
    source_material_id = '5'
    block = air
  []
  [Walls]
    type = FileNeutronicsMaterial
    transport_system = Neutron
    file_name = 'data/xs_macro.xml'
    source_material_id = '7'
    block = walls
  []
[]

[Executioner]
  type = Steady
  solve_type = PJFNK
  petsc_options_iname = '-ksp_gmres_restart'
  petsc_options_value = ' 100'

  l_max_its = 100
  nl_abs_tol = 1e-12

  automatic_scaling = true
[]

[Outputs]
  exodus = true
[]
