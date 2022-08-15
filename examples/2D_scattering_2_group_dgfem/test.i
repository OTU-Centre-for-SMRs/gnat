# A simple test case with a purely absorbing medium and a point source in the
# middle of the domain.

[GlobalParams]
  num_groups = 2
  max_anisotropy = 0
  anisotropy = 0
[]

[Mesh]
  [domain]
    type = CartesianMeshGenerator
    dim = 2
    dx = 10
    dy = 10
    ix = 100
    iy = 100
  []
[]

[NeutronActivationStudy]
  [TransportSystem]
    scheme = upwinding_dfem
    execution_type = steady
    output_angular_fluxes = false

    order = FIRST
    family = MONOMIAL

    n_azimuthal = 2
    n_polar = 2

    vacuum_boundaries = 'left right top bottom'

    point_source_locations = '5.0 5.0 0.0'
    point_source_intensities = '1000.0'
    point_source_groups = '1'
  []
[]

[Materials]
  [Domain]
    type = ConstantNeutronicsMaterial
    group_absorption = '0.1 1.0'
    group_scattering = '0.5 1.0 0.0 0.0'
    group_speeds = '220000.0 220000.0'
  []
[]

[Problem]
  type = FEProblem
[]

[Executioner]
  type = Steady
  solve_type = PJFNK
  nl_rel_tol = 1e-12
[]
