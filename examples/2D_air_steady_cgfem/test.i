# A simple test case with a purely absorbing medium and a point source in the
# middle of the domain.

[Mesh]
  [domain]
    type = CartesianMeshGenerator
    dim = 2
    dx = '3.5 3 3.5'
    dy = '3.5 3 3.5'
    ix = '35 31 35'
    iy = '35 31 35'
    subdomain_id = '
      1 1 1
      1 1 1
      1 1 1'
  []
[]

[NeutronActivationStudy]
  execution_type = steady
  num_groups = 2
  max_anisotropy = 1

  [TransportSystem]
    scheme = saaf_cfem
    output_angular_fluxes = false

    order = FIRST
    family = LAGRANGE

    n_azimuthal = 3
    n_polar = 3

    vacuum_boundaries = 'left right top bottom'

    point_source_locations = '5.0 5.0 0.0'
    point_source_intensities = '1000.0'
    point_source_groups = '1'
  []
[]

[Materials]
  [Domain]
    type = FileNeutronicsMaterial
    num_groups = 2
    file_name = './examples/2D_air_steady_cgfem/cross_sections/cross_sections.txt'
    source_material_id = '1'
    block = '1'
  []
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
