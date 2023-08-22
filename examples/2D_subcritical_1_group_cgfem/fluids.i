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
  exodus = true
[]
