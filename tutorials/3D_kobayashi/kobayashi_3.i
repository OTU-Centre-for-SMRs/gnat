[Mesh]
  [domain]
    type = CartesianMeshGenerator
    dim = 3
    dx = '10 20 10 20'
    dy = '10 40 10 40'
    dz = '10 20 10 20'
    ix = '1 2 1 2'
    iy = '1 4 1 4'
    iz = '1 2 1 2'
    subdomain_id = '
    2 1 1 1
    0 1 1 1
    0 0 0 1
    1 1 1 1

    1 1 1 1
    1 1 1 1
    1 1 0 1
    1 1 1 1

    1 1 1 1
    1 1 1 1
    1 1 0 1
    1 1 0 1

    1 1 1 1
    1 1 1 1
    1 1 1 1
    1 1 1 1'
  []
  uniform_refine = 1
[]

[TransportSystems]
  [Neutron]
    num_groups = 1
    max_anisotropy = 0
    scheme = saaf_cfem
    particle_type = neutron

    order = FIRST
    family = LAGRANGE

    n_azimuthal = 5
    n_polar = 5

    vacuum_boundaries = 'right top front'
    reflective_boundaries = 'back bottom left'

    volumetric_source_blocks = '2'
    volumetric_source_moments = '1.0'
    volumetric_source_anisotropies = '0'
  []
[]

[TransportMaterials]
  [Shield]
    type = ConstantTransportMaterial
    transport_system = 'Neutron'
    group_total = '0.1'
    group_scattering = '0.05'
    anisotropy = 0
    block = '1 2'
  []
  [Duct]
    type = ConstantTransportMaterial
    transport_system = 'Neutron'
    group_total = '1.0e-4'
    group_scattering = '0.5e-4'
    anisotropy = 0
    block = '0'
  []
[]

[Executioner]
  type = Steady
  solve_type = PJFNK
  petsc_options_iname = '-pc_type -pc_hypre_type -ksp_gmres_restart'
  petsc_options_value = ' hypre    boomeramg      100'
  l_max_its = 50
  nl_rel_tol = 1e-12

  automatic_scaling = true
  off_diagonals_in_auto_scaling = true
  compute_scaling_once = false
[]

[VectorPostprocessors]
  [Line_A]
    type = LineValueSampler
    variable = flux_moment_1_0_0
    start_point = '5.0 5.0  5.0'
    end_point   = '5.0 95.0 5.0'
    num_points = 10
    sort_by = y
  []

  [Line_B]
    type = LineValueSampler
    variable = flux_moment_1_0_0
    start_point = '5.0  55.0 5.0'
    end_point   = '55.0 55.0 5.0'
    num_points = 6
    sort_by = x
  []

  [Line_C]
    type = LineValueSampler
    variable = flux_moment_1_0_0
    start_point = '5.0  55.0 35.0'
    end_point   = '55.0 55.0 35.0'
    num_points = 6
    sort_by = x
  []
[]

[Outputs]
  exodus = true
  csv = true
  execute_on = TIMESTEP_END
[]
