# A simple test case with a purely absorbing medium and a point source in the
# middle of the domain.

[Mesh]
  [domain]
    type = CartesianMeshGenerator
    dim = 1
    dx = 10
    ix = 100
  []
[]

[NeutronActivationStudy]
  [TransportSystem]
    scheme = upwinding_dfem
    execution_type = transient
    num_groups = 1
    output_angular_fluxes = true

    order = FIRST
    family = MONOMIAL
    constant_ic = 0.0

    n_azimuthal = 1
    n_polar = 1

    max_anisotropy = 0
    vacuum_boundaries = 'left right'

    point_source_locations = '5.0 0.0 0.0'
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
    group_absorption = 1.0
    group_speeds = 220000.0
  []
[]

[Problem]
  type = FEProblem
[]

[Executioner]
  type = Transient
  num_steps = 15
  solve_type = NEWTON

  [./TimeIntegrator]
    type = ImplicitEuler
  []

  [./TimeStepper]
    type = ConstantDT
    dt = 0.000004545
  []
[]
