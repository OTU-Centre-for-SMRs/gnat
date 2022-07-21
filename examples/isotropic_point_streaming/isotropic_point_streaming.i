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

    num_groups = 1
    max_anisotropy = 0
    vacuum_boundaries = 'left right top bottom'
    debug_disable_scattering = true

    point_source_locations = '0.0 0.0 0.0'
    point_source_intensities = '1000.0'
    point_source_groups = '1'
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
  type = FEProblem # This is the "normal" type of Finite Element Problem in MOOSE
  coord_type = XYZ
[]

[Executioner]
  type = Steady
  solve_type = NEWTON
[]

[Outputs]
  execute_on = 'timestep_end'
  exodus = true
[]
