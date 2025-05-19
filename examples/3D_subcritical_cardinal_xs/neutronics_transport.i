[Mesh]
  [NeutronicsDomain]
    type = FileMeshGenerator
    file = openmc_out.e
    use_for_exodus_restart = true
  []
  uniform_refine = 0
[]

[TransportSystems]
  [Neutron]
    num_groups = 3
    max_anisotropy = 0
    scheme = saaf_cfem
    particle_type = neutron

    flux_moment_names = 'Flux_Moment'
    debug_disable_fission = false
    eigen = true
    constant_ic = '1.0'

    order = FIRST
    family = LAGRANGE

    n_azimuthal = 3
    n_polar = 3

    vacuum_boundaries = '10000'
    reflective_boundaries = '10001'

    #volumetric_source_blocks = '5'
    #volumetric_source_anisotropies = '0'
    #volumetric_source_moments = '1.0 0.0 0.0 0.0'
  []
[]

[CardinalMGXS]
  xs_source = 'subapp'
  from_multi_app = 'cardinal_xs'

  transport_system = 'Neutron'
  scatter_anisotropy = 0
  add_fission_heating = true

  block = 'fuel cladding graphite_core graphite_reflector stainless_steel
           fuel_trimmer_tri graphite_core_trimmer_tri graphite_reflector_trimmer_tri
           stainless_steel_trimmer_tri'
[]

[TransportMaterials]
  [Void]
    type = VoidTransportMaterial
    block = 'air src src_trimmer_tri'
  []
[]

[Postprocessors]
  [TotalFlux]
    type = TotalFluxPostProcessor
    num_groups = 3
    group_scalar_fluxes = 'Flux_Moment_1_0_0 Flux_Moment_2_0_0 Flux_Moment_3_0_0'
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

  initial_eigenvalue = 1.2
  normal_factor = 1.0
  normalization = 'TotalFlux'
  free_power_iterations = 6

  solve_type = PJFNK
  #petsc_options_iname = '-pc_type -pc_hypre_type -ksp_gmres_restart'
  #petsc_options_value = ' hypre    boomeramg      600'
  petsc_options_iname = '-ksp_gmres_restart'
  petsc_options_value = ' 600'

  nl_abs_tol = 1e-12

  automatic_scaling = true
  off_diagonals_in_auto_scaling = true
  compute_scaling_once = false
[]

[MultiApps]
  [cardinal_xs]
    type = FullSolveMultiApp
    app_type = 'CardinalApp'
    execute_on = INITIAL
    input_files = 'openmc.i'
  []
[]

[Outputs]
  csv = true
  exodus = true
  execute_on = TIMESTEP_END
[]
