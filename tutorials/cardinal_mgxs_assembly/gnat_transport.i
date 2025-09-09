[Mesh]
  [NeutronicsDomain]
    type = FileMeshGenerator
    file = openmc_mgxs_out.e
    use_for_exodus_restart = true
  []
[]

[TransportSystems]
  [Neutron]
    num_groups = 2
    max_anisotropy = 0
    scheme = saaf_cfem
    particle_type = neutron

    flux_moment_names = 'Flux_Moment'
    debug_disable_fission = false
    eigen = true
    constant_ic = '1.0'

    order = FIRST
    family = LAGRANGE

    n_azimuthal = 5
    n_polar = 5

    vacuum_boundaries = '10002'
    reflective_boundaries = '10000 10001'
  []
[]

[CardinalMGXS]
  xs_source = 'mesh'
  transport_system = 'Neutron'
  scatter_anisotropy = 0
  add_fission_heating = true
[]

[Postprocessors]
  [TotalFlux]
    type = TotalFluxPostProcessor
    num_groups = 2
    group_scalar_fluxes = 'Flux_Moment_1_0_0 Flux_Moment_2_0_0'
    execute_on = LINEAR
    outputs = 'csv'
  []
[]

[VectorPostprocessors]
  [EigenValues]
    type = Eigenvalues
    inverse_eigenvalue = true
  []
[]

[Executioner]
  type = Eigenvalue

  initial_eigenvalue = 1.5
  normal_factor = 1.0
  normalization = 'TotalFlux'
  free_power_iterations = 6

  solve_type = PJFNK
  petsc_options_iname = '-ksp_gmres_restart'
  petsc_options_value = ' 600'

  nl_abs_tol = 1e-12

  automatic_scaling = true
  off_diagonals_in_auto_scaling = true
  compute_scaling_once = false
[]

[Outputs]
  csv = true
  exodus = true
  execute_on = 'TIMESTEP_END'
[]
