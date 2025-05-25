[DepletionLibrary]
  depletion_file = '../../data/depl/chain_endfb71_pwr.xml'
  depletion_file_source = openmc
  show_warnings = false
  add_data_userobject = true
[]

# Mesh is in cm, need to convert all units to cm.
[Mesh]
  [Domain]
    type = FileMeshGenerator
    file = flow_steady.e
    use_for_exodus_restart = true
  []
[]

[Physics]
  [NavierStokes]
    [Flow/flow]
      compressibility = 'incompressible'

      mass_advection_interpolation = 'average'
      momentum_advection_interpolation = 'upwind'

      # Rho and nu for air at 20C and atmospheric pressure.
      density = 0.001276 # g cm^{-3}
      dynamic_viscosity = 0.0001722 # g cm^{-1} s^{-2}

      initial_velocity = '1e-12 1e-12 0'
      initial_pressure = 0.0

      inlet_boundaries = 'outflow stack'
      momentum_inlet_types = 'fixed-velocity fixed-velocity'
      momentum_inlet_functors = '-100.0 0.0; 0.0 100.0'

      outlet_boundaries = 'inflow'
      momentum_outlet_types = 'fixed-pressure-zero-gradient'
      pressure_functors = '10.0'

      wall_boundaries = 'earth building sky'
      momentum_wall_types = 'noslip noslip slip'

      block = air
    []
    [Turbulence/ml]
      turbulence_handling = 'mixing-length'
      coupled_flow_physics = flow
      mixing_length_walls = 'earth building stack'
      mixing_length_aux_execute_on = 'initial'
    []
  []
[]

[TracerDepletionSystem]
  scheme = fv
  using_moose_ns_fv = true
  fv_adv_interpolation = upwind
  fv_face_interpolation = skewness-corrected
  turbulence_handling = mixing-length
  schmidt_number = 1.4
  temperature = 300.0

  extra_nuclides = 'Cs137'
  extra_nuclide_number_densities = '0.0'

  block = air

  inlet_boundaries = 'stack'
  #Inlet concentration of roughly 1 MBq.
  inlet_number_densities = '1.4e15'

  debug_filter_nuclides = 'Cs137 Ba137_m1 Ba137'

  add_photon_sources = true
  # Just the 662 keV gamma photon.
  photon_group_boundaries = '663000.0 661000.0'
[]

[Executioner]
  type = Transient
  solve_type = NEWTON
  petsc_options_iname = '-pc_type -pc_factor_shift_type'
  petsc_options_value = ' lu       NONZERO'

  [TimeStepper]
    type = ConstantDT
    dt = 5.0
  []
  end_time = 60

  nl_abs_tol = 1e-8
  nl_max_its = 20
  l_max_its = 100
  line_search = 'none'

  automatic_scaling = true
  off_diagonals_in_auto_scaling = true
  compute_scaling_once = false
[]

[Outputs]
  exodus = true
[]
