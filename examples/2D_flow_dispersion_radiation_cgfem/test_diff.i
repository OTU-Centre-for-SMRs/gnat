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
    nx = 101
    ny = 21
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

[TransportSystems]
  [Neutron]
    scheme = saaf_cfem
    particle_type = neutron
    num_groups = 2
    max_anisotropy = 0

    output_angular_fluxes = false
    flux_moment_names = 'Flux_Moment'

    order = FIRST
    family = LAGRANGE

    n_azimuthal = 2
    n_polar = 2

    vacuum_boundaries = 'left right top bottom'

    point_source_locations = '5.0 0.0 0.0'
    point_source_intensities = '1000.0'
    point_source_groups = '1'

    constant_ic = '0.0 0.0'
  []
[]

[NuclideSystem]
  velocity_type = variable
  u_var = vel_x
  v_var = vel_y

  nuclides = 'H1 H2 H3'

  transport_system = Neutron

  order = FIRST
  family = LAGRANGE

  xs_file_name = './examples/2D_flow_dispersion_radiation_cgfem/xs_micro/cross_sections.txt'
  xs_type = micro
  xs_source_material_id = '1'

  nuclide_prop_file_name = './examples/2D_flow_dispersion_radiation_cgfem/nuclide_system_diff.txt'
  half_life_units = years

  density = 0.0

  [AddNuclideBCs]
    [H1_Inflow]
      type = ADIsotopeInflowBC
      boundary = left
      inflow_rate = 0.11186823

      excluded_nuclides = 'H2 H3'
    []

    [H2_Inflow]
      type = ADIsotopeInflowBC
      boundary = left
      inflow_rate = 0.000034823

      excluded_nuclides = 'H1 H3'
    []

    [Outflow]
      type = ADIsotopeOutflowBC
      boundary = right
    []
  []
[]

[TransportMaterials]
  [Water]
    type = FileNeutronicsMaterial
    transport_system = Neutron
    file_name = './examples/2D_flow_dispersion_radiation_cgfem/xs_macro/cross_sections.txt'
    source_material_id = '1'
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
