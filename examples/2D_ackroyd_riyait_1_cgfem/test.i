[Mesh]
  [domain]
    type = CartesianMeshGenerator
    dim = 2
    dx = '11 3 3 11'
    dy = '15 3 3 15'
    ix = '11 3 3 11'
    iy = '15 3 3 15'
    subdomain_id = '
    2 1 1 2
    2 0 0 2
    2 0 0 2
    2 1 1 2'
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

    n_azimuthal = 13
    n_polar = 13

    vacuum_boundaries = 'left right top bottom'

    volumetric_source_blocks = '0'
    volumetric_source_moments = '1.0'
    volumetric_source_anisotropies = '0'
    scale_sources = true
  []
[]

[TransportMaterials]
  [Duct]
    type = AbsorbingNeutronicsMaterial
    transport_system = Neutron
    group_absorption = '0.0'
    group_speeds = '220000'
    block = '1'
  []
  [Other]
    type = AbsorbingNeutronicsMaterial
    transport_system = Neutron
    group_absorption = '0.5'
    group_speeds = '220000'
    block = '0 2'
  []
[]

[Executioner]
  type = Steady
  solve_type = PJFNK
  #petsc_options_iname = '-pc_type -pc_hypre_type -ksp_gmres_restart'
  #petsc_options_value = ' hypre    boomeramg      10'
  l_max_its = 50
  nl_rel_tol = 1e-12

  automatic_scaling = true
  off_diagonals_in_auto_scaling = true
  compute_scaling_once = false
[]

[Outputs]
  exodus = true
[]
