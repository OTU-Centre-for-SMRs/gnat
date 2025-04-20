[Mesh]
  [NeutronicsDomain]
    type = FileMeshGenerator
    file = m_and_ms.e
  []
[]

[TransportSystems]
  [Neutron]
    num_groups = 7
    max_anisotropy = 0
    scheme = saaf_cfem
    particle_type = neutron

    order = FIRST
    family = LAGRANGE

    n_azimuthal = 3
    n_polar = 3

    flux_moment_names = 'Flux_Moment'
    debug_disable_fission = false
    eigen = true
    constant_ic = '1.0'

    reflective_boundaries = 'vacuum'
  []
[]

[TransportMaterials]
  [Letters]
    type = FileTransportMaterial
    transport_system = Neutron
    file_name = 'C5G7_XS.xml'
    source_material_id = '3' # MOX 8.7%
    block = 'letters'
  []
  [M_M]
    type = FileTransportMaterial
    transport_system = Neutron
    file_name = 'C5G7_XS.xml'
    source_material_id = '0' # UO2
    block = 'm1 m2'
  []
  [Other]
    type = FileTransportMaterial
    transport_system = Neutron
    file_name = 'C5G7_XS.xml'
    source_material_id = '6' # Moderator
    block = 'bckgrnd'
  []
[]

[Postprocessors]
  [TotalFlux]
    type = TotalFluxPostProcessor
    num_groups = 7
    group_scalar_fluxes = 'Flux_Moment_1_0_0 Flux_Moment_2_0_0 Flux_Moment_3_0_0 Flux_Moment_4_0_0 Flux_Moment_5_0_0 Flux_Moment_6_0_0 Flux_Moment_7_0_0'
    execute_on = LINEAR
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

  initial_eigenvalue = 1.0
  normal_factor = 1.0
  normalization = 'TotalFlux'
  free_power_iterations = 6

  solve_type = PJFNK
  #petsc_options_iname = '-pc_type -pc_hypre_type -ksp_gmres_restart'
  #petsc_options_value = ' hypre    boomeramg      600'
  petsc_options_iname = '-ksp_gmres_restart'
  petsc_options_value = ' 600'

  nl_abs_tol = 1e-8

  automatic_scaling = true
  off_diagonals_in_auto_scaling = true
  compute_scaling_once = false
[]

[Outputs]
  csv = true
  exodus = true
  execute_on = 'TIMESTEP_END'
[]
