# This input deck computes collided fluxes using the merged uncollided fluxes for a
# 2D version of the OTU subcritical assembly.

[Mesh]
  [NeutronicsDomain]
    type = FileMeshGenerator
    file = 2D_quarter_core.e
  []
  uniform_refine = 0
[]

[TransportSystems]
  [Neutron]
    num_groups = 4
    max_anisotropy = 0
    scheme = saaf_cfem
    particle_type = neutron
    flux_moment_names = 'Flux_Moment'
    debug_disable_fission = false

    order = FIRST
    family = LAGRANGE

    n_azimuthal = 3
    n_polar = 3

    vacuum_boundaries = 'vacuum'
    reflective_boundaries = 'reflective'

    uncollided_from_multi_app = 'Merge'
    use_conservative_uncollided_transfers = false
  []
[]

[TransportMaterials]
  [Air]
    type = FileTransportMaterial
    transport_system = Neutron
    file_name = 'macro_xs.xml'
    source_material_id = '5' #
    block = 'air_gap src_air'
  []

  [Cladding]
    type = FileTransportMaterial
    transport_system = Neutron
    file_name = 'macro_xs.xml'
    source_material_id = '7' #
    block = 'cladding'
  []

  [Fuel]
    type = FileTransportMaterial
    transport_system = Neutron
    file_name = 'macro_xs.xml'
    source_material_id = '6' #
    block = 'fuel'
  []

  [Graphite]
    type = FileTransportMaterial
    transport_system = Neutron
    file_name = 'macro_xs.xml'
    source_material_id = '8' #
    block = 'core_graphite reflector_graphite'
  []

  [Control]
    type = FileTransportMaterial
    transport_system = Neutron
    file_name = 'macro_xs.xml'
    source_material_id = '8' #
    block = 'control_or_graphite'
  []
[]

[Executioner]
  type = Steady
  solve_type = PJFNK
  petsc_options_iname = '-pc_type -pc_hypre_type -ksp_gmres_restart'
  petsc_options_value = ' hypre    boomeramg      600'

  l_max_its = 50
  nl_abs_tol = 1e-12

  automatic_scaling = true
  off_diagonals_in_auto_scaling = true
  compute_scaling_once = false
[]

[Outputs]
  exodus = true
  execute_on = 'TIMESTEP_END'
[]

[MultiApps]
  [Merge]
    type = FullSolveMultiApp
    app_type = GnatApp
    input_files = 'neutronics_merge.i'
    execute_on = INITIAL
  []
[]
