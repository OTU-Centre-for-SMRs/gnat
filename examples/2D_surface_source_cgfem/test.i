[Mesh]
  [gen]
    type = GeneratedMeshGenerator
    dim = 2
    xmin = 0
    xmax = 10
    ymin = 0
    ymax = 10
    nx = 51
    ny = 51
  []
[]

[TransportSystems]
  [Neutron]
    scheme = saaf_cfem
    particle_type = neutron
    num_groups = 1
    max_anisotropy = 0

    output_angular_fluxes = false
    flux_moment_names = 'Flux_Moment'

    order = FIRST
    family = LAGRANGE

    n_azimuthal = 6
    n_polar = 6

    vacuum_boundaries = 'bottom'
    source_boundaries = 'top left right'

    boundary_source_anisotropy = '0 0 0 '
    boundary_source_moments = '1000; 1000; 1000'
  []
[]

[TransportMaterials]
  [Water]
    type = FileNeutronicsMaterial
    transport_system = Neutron
    file_name = './examples/2D_surface_source_cgfem/xs_macro/cross_sections.txt'
    source_material_id = '1'
  []

  #[Void]
  #  type = VoidNeutronicsMaterial
  #  transport_system = Neutron
  #  group_speeds = '220000 220000'
  #[]
[]

[Problem]
  type = FEProblem
[]

[Executioner]
  type = Steady
  solve_type = PJFNK
  petsc_options_iname = '-pc_type -pc_hypre_type -ksp_gmres_restart'
  petsc_options_value = 'hypre boomeramg 10'
  l_max_its = 50
  nl_rel_tol = 1e-12
[]

#[Executioner]
#  type = Steady
#  solve_type = NEWTON
#  petsc_options_iname = '-pc_type -pc_factor_shift_type'
#  petsc_options_value = 'lu       NONZERO'
#
#  nl_abs_tol = 1e-9
#  nl_max_its = 150
#  line_search = 'none'
#
#  automatic_scaling = true
#  off_diagonals_in_auto_scaling = true
#  compute_scaling_once = false
#[]
