# The 2D Ackroyd-Riyait voided dog-legged duct benchmark problem.

[Mesh]
  [NeutronicsDomain]
    type = FileMeshGenerator
    file = kirk_name.e
  []
[]

[TransportSystems]
  [Neutron]
    num_groups = 7
    max_anisotropy = 0
    scheme = saaf_cfem
    particle_type = neutron

    debug_disable_fission = false
    eigen = true
    constant_ic = '1.0'

    order = FIRST
    family = LAGRANGE

    n_azimuthal = 3
    n_polar = 3

    vacuum_boundaries = 'vacuum'
  []
[]

[TransportMaterials]
  [I]
    type = FileTransportMaterial
    transport_system = Neutron
    file_name = 'C5G7_XS.xml'
    source_material_id = '0'
    block = 'src'
  []
  [Other]
    type = FileTransportMaterial
    transport_system = Neutron
    file_name = 'C5G7_XS.xml'
    source_material_id = '6'
    block = 'shield'
  []
[]

[Executioner]
  type = Eigenvalue

  initial_eigenvalue = 1.0
  normal_factor = 1.0
  normalization = 'TotalFlux'
  free_power_iterations = 6

  solve_type = PJFNK
  petsc_options_iname = '-ksp_gmres_restart'
  petsc_options_value = ' 600'

  nl_abs_tol = 1e-8

  automatic_scaling = true
  off_diagonals_in_auto_scaling = true
  compute_scaling_once = false
[]

[Postprocessors]
  [TotalFlux]
    type = TotalFluxPostProcessor
    num_groups = 7
    group_scalar_fluxes = 'flux_moment_1_0_0 flux_moment_2_0_0 flux_moment_3_0_0 flux_moment_4_0_0 flux_moment_5_0_0 flux_moment_6_0_0 flux_moment_7_0_0'
    execute_on = LINEAR
  []
[]

[VectorPostprocessors]
  [EigenValues]
    type = Eigenvalues
    inverse_eigenvalue = true
  []
[]

[Outputs]
  csv = true
  exodus = true
  execute_on = 'TIMESTEP_END'
[]
