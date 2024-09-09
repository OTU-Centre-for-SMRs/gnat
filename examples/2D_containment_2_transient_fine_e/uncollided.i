[Mesh]
  [NeutronicsDomain]
    type = FileMeshGenerator
    file = Neutronics_Mesh.e
  []
  uniform_refine = 2
[]

[UncollidedFlux]
  [Neutron]
    uncollided_flux_treatment = sasf
    num_groups = 8
    max_anisotropy = 0

    point_source_locations = '0.0 0.0 0.0'
    point_source_moments = '7.534300e0 1.060700e0 3.250000e-2 6.347148e2 3.566577e2 0.0 0.0 0.0'
    point_source_anisotropies = '0'

    sasf_near_source_boundary = 'barrel'
    sasf_vacuum_boundaries = 'vacuum'
    sasf_near_source_cross_sections = '0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0'

    is_conservative_transfer_src = false
  []
[]

[TransportMaterials]
  [Stainless]
    type = FileTransportMaterial
    transport_system = Neutron
    file_name = 'data/xs_macro.xml'
    source_material_id = '6'
    block = stainless
  []
  [Air]
    type = FileTransportMaterial
    transport_system = Neutron
    file_name = 'data/xs_macro.xml'
    source_material_id = '5'
    block = air
  []
  [Concrete]
    type = FileTransportMaterial
    transport_system = Neutron
    file_name = 'data/xs_macro.xml'
    source_material_id = '7'
    block = concrete
  []
  [Water]
    type = FileTransportMaterial
    transport_system = Neutron
    file_name = 'data/water_macro_xs.xml'
    source_material_id = '1'
    block = light_water
  []
[]

[Executioner]
  type = Steady
  solve_type = PJFNK
  petsc_options_iname = '-ksp_gmres_restart'
  petsc_options_value = ' 600'

  l_max_its = 50
  nl_abs_tol = 1e-12
  #nl_rel_tol = 1e-13

  automatic_scaling = true
  off_diagonals_in_auto_scaling = true
  compute_scaling_once = false
[]

[Outputs]
  exodus = true
[]
