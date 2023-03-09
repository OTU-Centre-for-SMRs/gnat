[Mesh]
  [Domain]
    type = FileMeshGenerator
    file = containment.e
  []
[]

[Modules]
  [NavierStokesFV]
    compressibility = 'incompressible'
    momentum_advection_interpolation = 'skewness-corrected'
    momentum_face_interpolation = 'skewness-corrected'
    pressure_face_interpolation = 'skewness-corrected'

    density = 1.276
    dynamic_viscosity = 1.0

    initial_velocity = '0 0 0'
    initial_pressure = 0.0

    inlet_boundaries = 'inlet'
    momentum_inlet_types = 'fixed-velocity'
    momentum_inlet_function = '0 -100.0' # cm/s (1 m/s)
    wall_boundaries = 'walls'
    momentum_wall_types = 'noslip'
    outlet_boundaries = 'outlet'
    momentum_outlet_types = 'fixed-pressure'
    pressure_function = '0'

    turbulence_handling = mixing-length
    mixing_length_walls = 'walls'

    block = air
  []
[]

[AuxVariables]
  [eddy_diff]
    type = MooseVariableFVReal
  []
[]

[AuxKernels]
  [EddyDiffusivity]
    type = MLEddyDiffusivity
    variable = eddy_diff
    u = vel_x
    v = vel_y
    mixing_length = mixing_length
  []
[]

[TransportSystems]
  [Neutron]
    scheme = flux_moment_transfer
    particle_type = neutron
    num_groups = 2

    from_multi_app = Neutronics
    from_flux_moment_names = Flux_Moment
    use_copy = false

    flux_moment_names = 'Flux_Moment'

    order = FIRST
    family = LAGRANGE
  []
[]

[NuclideSystem]
  velocity_type = variable
  u_var = vel_x
  v_var = vel_y
  eddy_diffusivity = eddy_diff

  nuclides = 'Ar40 Ar41 K41'

  transport_system = Neutron

  order = FIRST
  family = LAGRANGE

  xs_file_name = './examples/2D_containment/xs_micro/cross_sections.txt'
  xs_type = micro
  xs_source_material_id = '3'

  nuclide_prop_file_name = './examples/2D_containment/nuclide_system.txt'
  half_life_units = minutes

  block = air

  density = 0.0

  [AddNuclideBCs]
    [Ar40]
      type = ADIsotopeInflowBC
      boundary = inlet
      inflow_rate = 0.000011163

      excluded_nuclides = 'Ar41 K41'
    []

    [Outflow]
      type = ADIsotopeOutflowBC
      boundary = outlet
    []
  []
[]

[Executioner]
  type = Transient
  solve_type = 'NEWTON'
  petsc_options_iname = '-pc_type -pc_factor_shift_type'
  petsc_options_value = 'lu       NONZERO'

  [TimeStepper]
    type = ConstantDT
    dt = 0.125
  []
  end_time = 100.0

  nl_abs_tol = 1e-8
  nl_max_its = 50
  line_search = 'none'

  automatic_scaling = true
  off_diagonals_in_auto_scaling = true
  compute_scaling_once = false
[]

#[Outputs]
#  exodus = true
#[]

[MultiApps]
  [Neutronics]
    type = FullSolveMultiApp
    app_type = GnatApp
    input_files = 'neutronics.i'
    execute_on = INITIAL
  []
[]
