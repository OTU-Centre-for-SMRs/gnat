[Mesh]
  [NeutronicsDomain]
    type = FileMeshGenerator
    file = Radiation_Mesh.e
  []
[]

[TransportSystems]
  [Neutron]
    num_groups = 8
    max_anisotropy = 1
    scheme = saaf_cfem
    particle_type = neutron

    flux_moment_names = 'Flux_Moment'
    is_conservative_transfer_src = true

    order = FIRST
    family = LAGRANGE

    n_azimuthal = 4
    n_polar = 4

    vacuum_boundaries = 'vacuum'

    volumetric_source_blocks = fuel
    volumetric_source_moments = '
    40.5717 0.0 0.0 0.0
    0.0     0.0 0.0 0.0'
    volumetric_source_anisotropies = '0'
  []
[]

[TransportMaterials]
  [Air]
    type = FileTransportMaterial
    transport_system = Neutron
    file_name = 'data/macro_xs.xml'
    source_material_id = '5'
    block = air
  []
  [Cladding]
    type = FileTransportMaterial
    transport_system = Neutron
    file_name = 'data/macro_xs.xml'
    source_material_id = '6'
    block = cladding
  []
  [Fuel]
    type = FileTransportMaterial
    transport_system = Neutron
    file_name = 'data/macro_xs.xml'
    source_material_id = '7'
    block = fuel
  []
  [Wood]
    type = FileTransportMaterial
    transport_system = Neutron
    file_name = 'data/macro_xs.xml'
    source_material_id = '8'
    block = wood
  []
[]

[Executioner]
  type = Steady
  solve_type = PJFNK
  #petsc_options_iname = '-pc_type -pc_hypre_type -ksp_gmres_restart'
  #petsc_options_value = ' hypre    boomeramg      600'
  petsc_options_iname = '-ksp_gmres_restart'
  petsc_options_value = ' 600'

  l_max_its = 50
  nl_abs_tol = 1e-12

  automatic_scaling = true
  off_diagonals_in_auto_scaling = true
  compute_scaling_once = false
[]

[Outputs]
  exodus = true
  execute_on = 'TIMESTEP_END'
[]
