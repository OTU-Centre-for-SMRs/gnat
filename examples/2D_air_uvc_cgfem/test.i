# A simple test case with a purely absorbing medium and a point source in the
# middle of the domain.

[Mesh]
  [domain]
    type = CartesianMeshGenerator
    dim = 2
    dx = 300.0
    dy = 300.0
    ix = 101
    iy = 101
  []
[]

[TransportSystems]
  [UVC]
    num_groups = 1
    scheme = saaf_cfem
    particle_type = photon

    order = FIRST
    family = LAGRANGE

    n_azimuthal = 18
    n_polar = 18

    max_anisotropy = 0
    vacuum_boundaries = 'left right top bottom'

    point_source_locations = '298.0 298.0 0.0'
    point_source_moments = '0.027646015'
    point_source_anisotropies = '0'
    scale_sources = true
  []
[]

[TransportMaterials]
  [Domain]
    type = ConstantNeutronicsMaterial
    transport_system = UVC
    group_absorption = 2.83e-5
    group_scattering = 4.6e-6
    group_speeds = 0.0 # Need to find a way to ignore params if photon. Might need to handle it manually.
  []
[]

[Executioner]
  type = Steady
  solve_type = PJFNK
  petsc_options_iname = '-pc_type -pc_factor_shift_type'
  petsc_options_value = ' lu       NONZERO'

  l_max_its = 50
  nl_abs_tol = 1e-12

  automatic_scaling = true
  off_diagonals_in_auto_scaling = true
  compute_scaling_once = false
[]

[Outputs]
  exodus = true
[]
