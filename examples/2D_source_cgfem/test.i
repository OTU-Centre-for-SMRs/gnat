# A simple test case with a purely absorbing medium and a point source in the
# middle of the domain.

[Mesh]
  [domain]
    type = CartesianMeshGenerator
    dim = 2
    dx = '4 2 4'
    dy = '4 2 4'
    ix = '40 20 40'
    iy = '40 20 40'
    subdomain_id = '
      1 1 1
      1 2 1
      1 1 1'
  []
[]

[NeutronActivationStudy]
  [TransportSystem]
    scheme = saaf_cfem
    execution_type = steady
    num_groups = 1
    output_angular_fluxes = true

    order = FIRST
    family = LAGRANGE

    n_azimuthal = 2
    n_polar = 2

    max_anisotropy = 0
    vacuum_boundaries = 'left right top bottom'

    debug_verbosity = level0
  []
[]

[Materials]
  [Domain]
    type = ConstantNeutronicsMaterial
    num_groups = 1
    group_absorption = 1.0
    anisotropy = 0
    group_scattering = 1.0
    group_speeds = 220000.0
    block = '1'
  []
  [Source]
    type = SourceNeutronicsMaterial
    num_groups = 1
    group_absorption = 1.0
    source_anisotropy = 0
    group_source = 1
    anisotropy = 0
    group_scattering = 1.0
    group_speeds = 220000.0
    block = '2'
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
