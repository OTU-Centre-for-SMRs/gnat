[Mesh]
  [NeutronicsDomain]
    type = FileMeshGenerator
    file = Radiation_Mesh.e
  []
[]

[TransportSystems]
  [Neutron]
    num_groups = 4
    max_anisotropy = 1
    scheme = saaf_cfem
    particle_type = neutron

    flux_moment_names = 'Flux_Moment'
    debug_disable_fission = false
    eigen = true

    order = FIRST
    family = LAGRANGE

    n_azimuthal = 4
    n_polar = 4

    vacuum_boundaries = 'vacuum'
  []
[]

[TransportMaterials]
  [Air]
    type = FileNeutronicsMaterial
    transport_system = Neutron
    file_name = 'macro_xs.xml'
    source_material_id = '5' #
    block = '2 7 8 9 '
  []

  [Cladding]
    type = FileNeutronicsMaterial
    transport_system = Neutron
    file_name = 'macro_xs.xml'
    source_material_id = '7' #
    block = 4
  []

  [Fuel]
    type = FileNeutronicsMaterial
    transport_system = Neutron
    file_name = 'macro_xs.xml'
    source_material_id = '6' #
    block = 3
  []

  [Graphite]
    type = FileNeutronicsMaterial
    transport_system = Neutron
    file_name = 'macro_xs.xml'
    source_material_id = '8' #
    block = '10 11'
  []

  [Control]
    type = FileNeutronicsMaterial
    transport_system = Neutron
    file_name = 'macro_xs.xml'
    source_material_id = '9' #
    block = 6
  []
[]

[Outputs]
  csv = true
  exodus = true
  execute_on = 'TIMESTEP_END'
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
  free_power_iterations = 6

  solve_type = PJFNK
  petsc_options_iname = '-pc_type -pc_hypre_type -ksp_gmres_restart'
  petsc_options_value = ' hypre    boomeramg      300'

  nl_abs_tol = 1e-12

  automatic_scaling = true
  off_diagonals_in_auto_scaling = true
  compute_scaling_once = false
[]
