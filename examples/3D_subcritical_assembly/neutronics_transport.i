[GlobalParams]
  file_name = '../../data/mgxs/full_sc/rodded/group_3_L_0.xml'
[]

[Mesh]
  [NeutronicsDomain]
    type = FileMeshGenerator
    file = mesh_out.e
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

[TransportMaterials]
  [Air]
    type = FileTransportMaterial
    transport_system = Neutron
    source_material_id = '0'
    block = '2 6 106'
  []

  [Cladding]
    type = FileTransportMaterial
    transport_system = Neutron
    source_material_id = '1'
    block = '1'
  []

  [Control]
    type = FileTransportMaterial
    transport_system = Neutron
    source_material_id = '2'
    block = '5 105'
  []

  [Fuel]
    type = FileTransportMaterial
    transport_system = Neutron
    source_material_id = '3'
    block = '0 100'
    add_heating = true
  []

  [CoreGraphite]
    type = FileTransportMaterial
    transport_system = Neutron
    source_material_id = '4'
    block = '3 103'
  []

  [ReflectorGraphite]
    type = FileTransportMaterial
    transport_system = Neutron
    source_material_id = '5'
    block = '4 104'
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

[Outputs]
  csv = true
  exodus = true
  execute_on = 'TIMESTEP_END'
[]
