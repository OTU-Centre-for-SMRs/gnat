# A simple test case with a purely absorbing medium and a point source in the
# middle of the domain.

[Mesh]
  [domain]
    type = CartesianMeshGenerator
    dim = 2
    dx = '11 2.0 1.0 2.0 11'
    dy = '11 3.0 3.0 3.0 3.0 11'
    ix = '11 2 1 2 11'
    iy = '11 3 3 3 3 11'
    subdomain_id = '
    1 0 0 0 1
    1 0 0 0 1
    1 0 0 0 1
    1 0 0 0 1
    1 0 0 0 1
    1 0 0 0 1'
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

    n_azimuthal = 10
    n_polar = 10

    vacuum_boundaries = 'left right top bottom'

    point_source_groups = '1'
    point_source_intensities = '1'
    point_source_locations = '14.0 17.0 0.0'

    debug_disable_scattering = true
  []
[]

[TransportMaterials]
  [Shield]
    type = AbsorbingNeutronicsMaterial
    transport_system = ''
    group_absorption = '1.0'
    group_speeds = '220000'
    block = '1'
  []
  [Duct]
    type = AbsorbingNeutronicsMaterial
    transport_system = ''
    group_absorption = '0.0'
    group_speeds = '220000'
    block = '0'
  []
  #[Other]
  #  type = AbsorbingNeutronicsMaterial
  #  transport_system = ''
  #  group_absorption = '0.5'
  #  group_speeds = '220000'
  #  block = '0'
  #[]
[]

[Problem]
  type = FEProblem
[]

[Outputs]
  exodus = true
[]

[Executioner]
  type = Steady
  solve_type = PJFNK
  #petsc_options_iname = '-pc_type -pc_hypre_type -ksp_gmres_restart'
  #petsc_options_value = 'hypre boomeramg 10'
  #petsc_options_iname = '-pc_type -pc_factor_shift_type'
  #petsc_options_value = 'lu       NONZERO'

  l_max_its = 50
  nl_abs_tol = 1e-12

  automatic_scaling = true
  off_diagonals_in_auto_scaling = true
  compute_scaling_once = false
[]
