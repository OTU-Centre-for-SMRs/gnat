[Mesh]
  [domain]
    type = FileMeshGenerator
    file = 'modified_kobayashi_3.e'
  []
  uniform_refine = 1
[]

[UncollidedFlux]
  [Neutron]
    uncollided_flux_treatment = sasf
    num_groups = 1
    max_anisotropy = 1

    point_source_locations = '0.0 0.0 0.0'
    point_source_moments = '1.0'
    point_source_anisotropies = '0'

    sasf_near_source_boundary = 'analytical_flux'
    sasf_vacuum_boundaries = 'vacuum'
    sasf_near_source_cross_sections = '0.1'
  []
[]

[Materials]
  [Duct]
    type = AbsorbingTransportMaterial
    transport_system = 'Neutron'
    num_groups = 1
    group_total = '1e-4'
    saaf_eta = 1.0
    block = duct
  []
  [Shield]
    type = AbsorbingTransportMaterial
    transport_system = 'Neutron'
    num_groups = 1
    group_total = '0.1'
    saaf_eta = 1.0
    block = shield
  []
[]

[Executioner]
  type = Steady
  solve_type = PJFNK
  petsc_options_iname = '-pc_type -pc_hypre_type -ksp_gmres_restart'
  petsc_options_value = ' hypre    boomeramg      30'

  l_max_its = 50
  nl_abs_tol = 1e-12

  automatic_scaling = true
  off_diagonals_in_auto_scaling = true
  compute_scaling_once = false
[]

[Outputs]
  exodus = true
  execute_on = TIMESTEP_END
[]

[Postprocessors]
  [H]
    type = AverageElementMinSize
  []
[]
