# A simple test case with a purely absorbing medium and a point source in the
# middle of the domain.

[Mesh]
  [domain]
    type = CartesianMeshGenerator
    dim = 3
    dx = 10
    dy = 10
    dz = 10
    ix = 20
    iy = 20
    iz = 20
  []
[]

[NeutronActivationStudy]
  [TransportSystem]
    scheme = upwinding_dfem
    execution_type = steady
    num_groups = 1
    output_angular_fluxes = true

    order = FIRST
    family = MONOMIAL

    n_azimuthal = 2
    n_polar = 2

    max_anisotropy = 0
    vacuum_boundaries = 'left right top bottom'

    point_source_locations = '5.0 5.0 5.0'
    point_source_intensities = '1000.0'
    point_source_groups = '1'

    debug_verbosity = level0
    debug_disable_scattering = true
  []
[]

[Materials]
  [Domain]
    type = AbsorbingNeutronicsMaterial
    num_groups = 1
    group_absorption = 0.0
    group_speeds = 2200.0
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
