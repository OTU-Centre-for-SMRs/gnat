# isotropic_point_streaming.i: A simple test case with a purely absorbing
# medium and a point source at the origin.

[Mesh]
  [domain]
    type = CartesianMeshGenerator
    dim = 2
    dx = 10
    dy = 10
    ix = 10
    iy = 10
  []
[]

[NeutronActivationStudy]
  [TransportSystem]
    execution_type = steady
    family = MONOMIAL
    order = FIRST

    n_azimuthal = 3
    n_polar = 3
    num_groups = 1

    max_anisotropy = 0
    vacuum_boundaries = 'left right top bottom'

    point_source_locations = '0.5 0.5 0.0'
    point_source_intensities = '1000.0'
    point_source_groups = '1'

    debug_disable_scattering = true
    debug_steady_state_ic = 1000.0
  []
[]

[Materials]
  [Domain]
    type = AbsorbingNeutronicsMaterial
    num_groups = 1
    group_removal = 1.0
    group_speeds = 1.0
  []
[]

[Problem]
  type = FEProblem
  coord_type = XYZ
[]

[Executioner]
  type = Steady
  solve_type = PJFNK
  petsc_options_iname = '-pc_type -pc_hypre_type'
  petsc_options_value = 'hypre boomeramg'
[]

[Outputs]
  exodus = true
[]
