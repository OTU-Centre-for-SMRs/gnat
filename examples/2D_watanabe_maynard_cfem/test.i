[Mesh]
  [domain]
    type = CartesianMeshGenerator
    dim = 2
    dx = '5.0 3.75 2.5 3.75 5.0'
    dy = '5.0 3.75 2.5 3.75 5.0'
    ix = '8 6 4 6 8'
    iy = '8 6 4 6 8'
    subdomain_id = '
    3 3 3 3 3
    3 2 2 2 3
    3 2 1 2 3
    3 2 2 2 3
    3 3 3 3 3'
  []
  uniform_refine = 2
[]

[TransportSystems]
  [Neutron]
    num_groups = 1
    max_anisotropy = 0
    scheme = saaf_cfem
    particle_type = neutron

    order = FIRST
    family = LAGRANGE

    n_azimuthal = 7
    n_polar = 7

    vacuum_boundaries = 'left right top bottom'

    volumetric_source_blocks = '1'
    volumetric_source_moments = '6.4'
    volumetric_source_anisotropies = '0'
    scale_sources = true
  []
[]

[TransportMaterials]
  [Void]
    type = ConstantTransportMaterial
    transport_system = Neutron
    group_total = '2.0e-10'
    group_scattering = '1.9e-10'
    block = '2'
  []
  [Other]
    type = ConstantTransportMaterial
    transport_system = Neutron
    group_total = '0.2'
    group_scattering = '0.19'
    block = '1 3'
  []
[]

[Executioner]
  type = Steady
  solve_type = PJFNK
  petsc_options_iname = '-ksp_gmres_restart'
  petsc_options_value = ' 300'

  l_max_its = 50
  nl_rel_tol = 1e-12

  automatic_scaling = true
  off_diagonals_in_auto_scaling = true
  compute_scaling_once = false
[]

[Outputs]
  exodus = true
  execute_on = 'TIMESTEP_END'
[]
