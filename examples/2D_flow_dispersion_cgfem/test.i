mu = 1.1
rho = 1.1

[Mesh]
  [gen]
    type = GeneratedMeshGenerator
    dim = 2
    xmin = 0
    xmax = 10
    ymin = -1
    ymax = 1
    nx = 100
    ny = 20
  []
[]

[Modules]
  [NavierStokesFV]
    compressibility = 'incompressible'

    density = ${rho}
    dynamic_viscosity = ${mu}

    initial_velocity = '0 0 0'
    initial_pressure = 0.0

    inlet_boundaries = 'left'
    momentum_inlet_types = 'fixed-velocity'
    momentum_inlet_function = '1 0'
    wall_boundaries = 'top bottom'
    momentum_wall_types = 'noslip noslip'
    outlet_boundaries = 'right'
    momentum_outlet_types = 'fixed-pressure'
    pressure_function = '0'
  []
[]

[NuclideSystem]
  velocity_type = variable
  u_var = vel_x
  v_var = vel_y

  nuclides = 'Ar40'

  order = FIRST
  family = LAGRANGE

  xs_file_name = './examples/2D_mass_transport_system_steady_cgfem/xs_micro/cross_sections.txt'
  xs_type = micro
  xs_source_material_id = '1'

  nuclide_prop_file_name = './examples/2D_mass_transport_system_steady_cgfem/nuclide_system.txt'
  half_life_units = minutes

  density = 0.0

  [AddNuclideBCs]
    [Inflow]
      type = ADIsotopeInflowBC
      boundary = left
      inflow_rate = 1.0
    []

    [Outflow]
      type = ADIsotopeOutflowBC
      boundary = right
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
    dt = 1.0
  []
  end_time = 15

  nl_abs_tol = 1e-9
  nl_max_its = 50
  line_search = 'none'

  automatic_scaling = true
  off_diagonals_in_auto_scaling = true
  compute_scaling_once = false
[]

[Outputs]
  exodus = true
[]
