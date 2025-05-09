[Mesh]
  [domain]
    type = CartesianMeshGenerator
    dim = 1
    dx = '2.0 1.0 2.0 1.0 2.0'
    ix = '4 2 4 2 4'
    subdomain_id = '0 1 2 3 4'
  []
[]

[TransportSystems]
  [Neutron]
    num_groups = 1
    scheme = saaf_cfem
    particle_type = neutron

    order = FIRST
    family = LAGRANGE

    n_polar = 10

    vacuum_boundaries = 'right'
    reflective_boundaries = 'left'

    volumetric_source_blocks = '0 3'
    volumetric_source_moments = '50.0; 1.0'
    volumetric_source_anisotropies = '0 0'
  []
[]

[TransportMaterials]
  [One]
    type = ConstantTransportMaterial
    transport_system = Neutron
    group_total = '50.0'
    group_scattering = '0.0'
    block = '0'
  []
  [Two]
    type = ConstantTransportMaterial
    transport_system = Neutron
    group_total = '5.0'
    group_scattering = '0.0'
    block = '1'
  []
  [Three]
    type = ConstantTransportMaterial
    transport_system = Neutron
    group_total = '0.0'
    group_scattering = '0.0'
    block = '2'
  []
  [Four]
    type = ConstantTransportMaterial
    transport_system = Neutron
    group_total = '1.0'
    group_scattering = '0.9'
    block = '3'
  []
  [Five]
    type = ConstantTransportMaterial
    transport_system = Neutron
    group_total = '1.0'
    group_scattering = '0.9'
    block = '4'
  []
[]

[Executioner]
  type = Steady
  solve_type = PJFNK
  petsc_options_iname = '-ksp_gmres_restart'
  petsc_options_value = ' 100'
  l_max_its = 50
  nl_rel_tol = 1e-8

  automatic_scaling = true
  off_diagonals_in_auto_scaling = true
  compute_scaling_once = false
[]

[VectorPostprocessors]
  [scalar_flux]
    type = LineValueSampler
    sort_by = x
    variable = 'flux_moment_1_0_0'
    start_point = '0 0 0'
    end_point = '8 0 0'
    num_points = 1000
  []
[]

[Outputs]
  exodus = true
  csv = true
  execute_on = 'TIMESTEP_END'
[]
