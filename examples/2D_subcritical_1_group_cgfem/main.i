# Mesh is in cm, need to convert all units to cm.

[DepletionLibrary]
  depletion_file = 'data/chain_endfb71_pwr_air.xml'
  cross_section_file = 'data/air_micro_xs.xml'
  depletion_file_source = openmc
  show_warnings = false
[]

[Mesh]
  [Domain]
    type = FileMeshGenerator
    file = Flow_Mesh.e
  []
[]

[Modules]
  [NavierStokesFV]
    compressibility = 'incompressible'
    momentum_advection_interpolation = 'upwind'

    # Rho and nu for air at 20C and atmospheric pressure.
    density = 0.001276 # g cm^{-3}
    initial_temperature = 300.0
    dynamic_viscosity = 0.0001722 # g cm^{-1} s^{-2}

    inlet_boundaries = 'front_door back_door'
    momentum_inlet_types = 'fixed-velocity fixed-velocity'
    momentum_inlet_function = '-10.0 0.0; 10.0 0.0'

    wall_boundaries = 'walls'
    momentum_wall_types = 'noslip'

    outlet_boundaries = 'hvac'
    momentum_outlet_types = 'fixed-pressure-zero-gradient'
    pressure_function = '0.0'

    turbulence_handling = mixing-length
    mixing_length_walls = 'walls'

    block = air
  []
[]

[TransportSystems]
  [Neutron]
    scheme = flux_moment_transfer
    particle_type = neutron
    num_groups = 1

    from_multi_app = Neutronics
    from_flux_moment_names = Flux_Moment
    transfer_to_fv = true
    use_conservative_transfers = true

    flux_moment_names = 'Flux_Moment'

    order = CONSTANT
    family = MONOMIAL
  []
[]

[MobileDepletionSystem]
  scheme = fv
  using_moose_ns_fv = true
  fv_adv_interpolation = upwind
  fv_face_interpolation = skewness-corrected
  turbulence_handling = mixing-length
  schmidt_number = 1.4

  elements = 'N O Ar C'
  element_atom_fractions = '0.78084 0.210294 0.00934 0.000417'

  transport_system = Neutron

  block = air

  inlet_boundaries = 'front_door back_door'
  # Same inlet boundary element atom fraction as the initial mixture for testing purposes.
  inlet_atom_fractions = '0.78084 0.210294 0.00934 0.000417; 0.78084 0.210294 0.00934 0.000417'

  debug_filter_nuclides = 'Ar40 Ar41 K41'
[]

[AuxVariables]
  [Ar41Activity]
    order = CONSTANT
    family = MONOMIAL
    fv = true
  []
[]

[AuxKernels]
  [CompAr41Activity]
    type = SpecificActivity
    variable = Ar41Activity
    nuclide_fun = Ar41
    decay_const = 0.000105396
    is_fe = false
  []
[]

[Postprocessors]
  [Ar41ActOutletFlowrate]
    type = VolumetricFlowRateQp
    boundary = hvac
    u = vel_x
    v = vel_y

    advected_quantity = Ar41Activity
  []

  [Ar41ActTotal]
    type = ElementIntegralVariablePostprocessor
    variable = Ar41Activity
    block = air
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
  end_time = 3600.0
  steady_state_detection = true

  nl_abs_tol = 1e-8
  nl_max_its = 20
  l_max_its = 100
  line_search = 'none'

  automatic_scaling = true
  off_diagonals_in_auto_scaling = true
  compute_scaling_once = false
[]

[Outputs]
  csv = true
  exodus = true
[]

[MultiApps]
  [Neutronics]
    type = FullSolveMultiApp
    app_type = GnatApp
    input_files = 'neutronics.i'
    execute_on = INITIAL
  []
[]
