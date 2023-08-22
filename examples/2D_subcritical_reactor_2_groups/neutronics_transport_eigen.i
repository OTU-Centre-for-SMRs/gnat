[Mesh]
  [NeutronicsDomain]
    type = FileMeshGenerator
    file = Radiation_Mesh.e
  []
[]

[TransportSystems]
  [Neutron]
    num_groups = 2
    max_anisotropy = 0
    scheme = saaf_cfem
    particle_type = neutron

    output_angular_fluxes = false
    flux_moment_names = 'Flux_Moment'
    debug_disable_fission = false
    #use_scattering_jacobians = true
    eigen = true

    order = FIRST
    family = LAGRANGE

    n_azimuthal = 3
    n_polar = 3

    vacuum_boundaries = 'vacuum'
  []
[]

[TransportMaterials]
  [Air]
    type = FileNeutronicsMaterial
    transport_system = Neutron
    file_name = 'xs_macro/air_cross_sections.txt'
    source_material_id = '5' #
    block = '2 7 8 9 '
  []

  [Cladding]
    type = FileNeutronicsMaterial
    transport_system = Neutron
    file_name = 'xs_macro/clad_cross_sections.txt'
    source_material_id = '7' #
    block = 4
  []

  [Fuel]
    type = FileNeutronicsMaterial
    transport_system = Neutron
    file_name = 'xs_macro/fuel_cross_sections.txt'
    source_material_id = '6' #
    block = 3
  []

  [Graphite]
    type = FileNeutronicsMaterial
    transport_system = Neutron
    file_name = 'xs_macro/mod_cross_sections.txt'
    source_material_id = '8' #
    block = '10 11'
  []

  [Control]
    type = FileNeutronicsMaterial
    transport_system = Neutron
    file_name = 'xs_macro/steel_cross_sections.txt'
    source_material_id = '9' #
    block = 6
  []
[]

[Outputs]
  exodus = true
[]

[Executioner]
  type = Eigenvalue

  initial_eigenvalue = 1.2

  solve_type = PJFNK
  petsc_options_iname = '-pc_type -pc_factor_shift_type'
  petsc_options_value = 'lu       NONZERO'

  nl_abs_tol = 1e-12

  automatic_scaling = true
  off_diagonals_in_auto_scaling = true
  compute_scaling_once = false
[]
