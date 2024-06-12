[DepletionLibrary]
  depletion_file = 'data/chain_endfb71_pwr_air.xml'
  cross_section_file = 'data/air_xs_micro.xml'
  depletion_file_source = openmc
  show_warnings = false
[]

[Mesh]
  # Mesh is in cm, need to convert all units to cm.
  [Domain]
    type = FileMeshGenerator
    file = Flow_Mesh.e
    #use_for_exodus_restart = true
  []
[]

[Modules]
  [NavierStokesFV]
    compressibility = 'incompressible'
    momentum_advection_interpolation = 'upwind'

    # Rho and nu for air at 20C and atmospheric pressure.
    density = 0.001276 # g cm^{-3}
    dynamic_viscosity = 0.0001722 # g cm^{-1} s^{-2}

    initial_velocity = '0 0 0'
    initial_pressure = 0.0

    inlet_boundaries = 'inlet'
    momentum_inlet_types = 'fixed-velocity'
    momentum_inlet_function = '0 -100.0' # cm / s
    wall_boundaries = 'noslip'
    momentum_wall_types = 'noslip'
    outlet_boundaries = 'outlet'
    momentum_outlet_types = 'fixed-pressure'
    pressure_function = '0'

    turbulence_handling = mixing-length
    mixing_length_walls = 'noslip'

    block = air
  []
[]

[MobileDepletionSystem]
  scheme = supg_fe
  transport_system = Neutron

  order = FIRST
  family = LAGRANGE

  using_moose_ns_fv = true
  turbulence_handling = mixing-length
  schmidt_number = 1.4

  elements = 'N O Ar C'
  element_atom_fractions = '0.78084 0.210294 0.00934 0.000417'

  block = air

  outlet_boundaries = 'outlet'
  inlet_boundaries = 'inlet'
  inlet_atom_fractions = '0.78084 0.210294 0.00934 0.000417'

  debug_filter_nuclides = 'Ar40 Ar41 K41'
[]

[TransportSystems]
  [Neutron]
    scheme = flux_moment_transfer
    particle_type = neutron
    num_groups = 8

    from_multi_app = Neutronics
    from_flux_moment_names = Flux_Moment
    use_copy = false
    #init_from_file = true

    flux_moment_names = 'Flux_Moment'

    order = FIRST
    family = LAGRANGE
  []
[]

[AuxVariables]
  [Ar41Activity]
    order = FIRST
    family = LAGRANGE
  []
[]

[AuxKernels]
  [CompAr41Activity]
    type = SpecificActivity
    variable = Ar41Activity
    is_fe = true
    nuclide_var = Ar41
    decay_const = 0.000105396
  []
[]

[Executioner]
  type = Transient
  solve_type = NEWTON
  petsc_options_iname = '-pc_type -pc_factor_shift_type'
  petsc_options_value = 'lu       NONZERO'

  [TimeStepper]
    type = ConstantDT
    dt = 0.5
  []
  end_time = 100.0

  nl_abs_tol = 1e-8
  nl_max_its = 50
  line_search = 'none'

  automatic_scaling = true
  off_diagonals_in_auto_scaling = true
  compute_scaling_once = false
[]

[MultiApps]
  [Neutronics]
    type = FullSolveMultiApp
    app_type = GnatApp
    input_files = 'neutronics.i'
    execute_on = INITIAL
  []
[]

[Outputs]
  exodus = true
[]
