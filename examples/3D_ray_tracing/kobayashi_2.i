[Mesh]
  [domain]
    type = CartesianMeshGenerator
    dim = 3
    dx = '50 20 50'
    dy = '90 20 90'
    dz = '50 20 50'
    ix = '5 2 5'
    iy = '9 2 9'
    iz = '5 2 5'
    subdomain_id = '
    1 1 1
    1 1 1
    1 1 1

    1 0 1
    1 2 1
    1 0 1

    1 1 1
    1 1 1
    1 1 1'
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

    vacuum_boundaries = 'left right top bottom'

    uncollided_from_multi_app = 'Uncollided'
  []
[]

[TransportMaterials]
  [Shield]
    type = ConstantNeutronicsMaterial
    transport_system = 'Neutron'
    group_absorption = '0.05'
    group_scattering = '0.05'
    anisotropy = 0
    group_speeds = '220000'
    block = '1 2'
  []
  [Duct]
    type = ConstantNeutronicsMaterial
    transport_system = 'Neutron'
    group_absorption = '0.5e-4'
    group_scattering = '0.5e-4'
    anisotropy = 0
    group_speeds = '220000'
    block = '0'
  []
  #[Shield]
  #  type = AbsorbingNeutronicsMaterial
  #  transport_system = 'Neutron'
  #  group_absorption = '0.1'
  #  group_speeds = '220000'
  #  block = '1 2'
  #[]
  #[Duct]
  #  type = AbsorbingNeutronicsMaterial
  #  transport_system = 'Neutron'
  #  group_absorption = '1e-4'
  #  group_speeds = '220000'
  #  block = '0'
  #[]
[]

[Executioner]
  type = Steady
  solve_type = PJFNK
  #petsc_options_iname = '-pc_type -pc_hypre_type -ksp_gmres_restart'
  #petsc_options_value = ' hypre    boomeramg      100'
  l_max_its = 50
  nl_rel_tol = 1e-12

  automatic_scaling = true
  off_diagonals_in_auto_scaling = true
  compute_scaling_once = false
[]

[Outputs]
  exodus = true
[]

[MultiApps]
  [Uncollided]
    type = FullSolveMultiApp
    app_type = GnatApp
    input_files = 'kobayashi_2_rt.i'
    execute_on = INITIAL
  []
[]
