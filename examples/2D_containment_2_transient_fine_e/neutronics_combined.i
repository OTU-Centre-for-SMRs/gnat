[Mesh]
  [NeutronicsDomain]
    type = FileMeshGenerator
    file = Neutronics_Mesh.e
  []
  uniform_refine = 0
[]

[TransportSystems]
  [Neutron]
    num_groups = 8
    max_anisotropy = 0
    scheme = saaf_cfem
    particle_type = neutron

    flux_moment_names = 'Flux_Moment'
    scale_sources = true

    order = FIRST
    family = LAGRANGE

    is_conservative_transfer_src = true

    n_azimuthal = 5
    n_polar = 5

    vacuum_boundaries = 'vacuum'
    current_boundaries = 'barrel'

    uncollided_from_multi_app = 'Uncollided'
    use_conservative_uncollided_transfers = false
  []
[]

[TransportMaterials]
  [Stainless]
    type = FileTransportMaterial
    transport_system = Neutron
    file_name = '../../data/mgxs/containment_xs_macro.xml'
    source_material_id = '6'
    block = stainless
  []
  [Air]
    type = FileTransportMaterial
    transport_system = Neutron
    file_name = '../../data/mgxs/containment_xs_macro.xml'
    source_material_id = '5'
    block = air
  []
  [Concrete]
    type = FileTransportMaterial
    transport_system = Neutron
    file_name = '../../data/mgxs/containment_xs_macro.xml'
    source_material_id = '7'
    block = concrete
  []
  [Water]
    type = FileTransportMaterial
    transport_system = Neutron
    file_name = '../../data/mgxs/containment_xs_macro.xml'
    source_material_id = '1'
    block = light_water
  []
[]

[Executioner]
  type = Steady
  solve_type = PJFNK
  petsc_options_iname = '-pc_type -pc_hypre_type -ksp_gmres_restart'
  petsc_options_value = ' hypre    boomeramg      10'

  l_max_its = 100
  nl_abs_tol = 1e-12

  automatic_scaling = true
[]

[MultiApps]
  [Uncollided]
    type = FullSolveMultiApp
    input_files = 'uncollided.i'
    execute_on = INITIAL
  []
[]

[Outputs]
  exodus = true
  execute_on = 'TIMESTEP_END'
[]
