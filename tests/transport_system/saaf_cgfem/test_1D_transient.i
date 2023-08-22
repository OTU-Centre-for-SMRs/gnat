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

[TransportSystems]
  [Neutron]
    scheme = saaf_cfem
    particle_type = neutron
    num_groups = 1
    output_angular_fluxes = true

    order = FIRST
    family = LAGRANGE
    constant_ic = 0.0

    n_azimuthal = 1
    n_polar = 1

    max_anisotropy = 0
    vacuum_boundaries = 'left right'

    point_source_locations = '5.0 0.0 0.0'
    point_source_moments = '1000.0'
    point_source_anisotropies = '1'

    debug_verbosity = level0
    debug_disable_scattering = true
  []
[]

[TransportMaterials]
  [Domain]
    type = AbsorbingNeutronicsMaterial
    transport_system = Neutron
    group_absorption = 1.0
    group_speeds = 220000.0
  []
[]

[Problem]
  type = FEProblem
[]

[Outputs]
  exodus = true
[]

[Executioner]
  type = Transient
  num_steps = 15
  solve_type = NEWTON

  [TimeIntegrator]
    type = ImplicitEuler
  []

  [TimeStepper]
    type = ConstantDT
    dt = 0.000004545
  []
[]
