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
  num_groups = 1
  execution_type = steady
  debug_verbosity = level0

  [TransportSystem]
    scheme = saaf_cfem
    output_angular_fluxes = true

    order = FIRST
    family = LAGRANGE

    n_azimuthal = 1
    n_polar = 1

    vacuum_boundaries = 'left right'

    point_source_locations = '5.0 0.0 0.0'
    point_source_intensities = '1000.0'
    point_source_groups = '1'
  []

  [AddMobileIsotope]
    order = FIRST
    family = LAGRANGE

    isotope_name = "cs_137"

    diffusion_coefficient_base = 1.0
    half_life = 1.0
    half_life_units = minutes

    absorption_cross_sections = '1.0'

    velocity_type = constant
    constant_velocity = '1.0 0.0 0.0'
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
  solve_type = NEWTON
[]
