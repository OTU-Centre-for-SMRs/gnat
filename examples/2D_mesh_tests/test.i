[Mesh]
  [NeutronicsDomain]
    type = FileMeshGenerator
    file = Neutronics_Full_2.e
  []
[]

[TransportSystems]
  [Neutron]
    num_groups = 2
    max_anisotropy = 0
    scheme = saaf_cfem
    particle_type = neutron
    use_scattering_jacobians = false

    order = FIRST
    family = LAGRANGE

    n_azimuthal = 2
    n_polar = 2

    vacuum_boundaries = 'vacuum'
    source_boundaries = 'neutron'

    boundary_source_anisotropy = '0'
    boundary_source_moments = '1e12 1e10'
  []
[]

[TransportMaterials]
  [RPV]
    type = FileNeutronicsMaterial
    transport_system = Neutron
    file_name = './examples/2D_mesh_tests/xs_macro/rpv_cross_sections.txt'
    source_material_id = '2'
    block = rpv
  []
  [Air]
    type = FileNeutronicsMaterial
    transport_system = Neutron
    file_name = './examples/2D_mesh_tests/xs_macro/air_cross_sections.txt'
    source_material_id = '3'
    block = air
  []
  [Walls]
    type = FileNeutronicsMaterial
    transport_system = Neutron
    file_name = './examples/2D_mesh_tests/xs_macro/walls_cross_sections.txt'
    source_material_id = '4'
    block = walls
  []
[]

[Problem]
  type = FEProblem
[]

[Executioner]
  type = Steady
  solve_type = NEWTON
  petsc_options_iname = '-pc_type -pc_hypre_type -ksp_gmres_restart'
  petsc_options_value = 'hypre boomeramg 10'

  nl_abs_tol = 1e-12
  nl_max_its = 1000
  line_search = 'none'

  automatic_scaling = true
  off_diagonals_in_auto_scaling = true
  compute_scaling_once = false
[]
