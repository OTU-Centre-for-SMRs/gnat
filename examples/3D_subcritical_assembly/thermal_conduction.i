[Mesh]
  [NeutronicsDomain]
    type = FileMeshGenerator
    file = neutronics_transport_4_groups_out.e
    use_for_exodus_restart = true
  []
  uniform_refine = 0
[]

[TransportSystems]
  [Neutron]
    num_groups = 4
    max_anisotropy = 0
    scheme = flux_moment_transfer
    particle_type = neutron

    flux_moment_names = 'Flux_Moment'
    debug_disable_fission = false
    init_from_file = true

    order = FIRST
    family = LAGRANGE
  []
[]

[Variables]
  [T]
    order = FIRST
    family = LAGRANGE
    initial_condition = 300.0
  []
[]

[Kernels]
  [Diffusion]
    type = HeatConduction
    variable = T
  []

  [Fission_Heating]
    type = FissionHeatSource
    variable = T
    num_groups = 4
    transport_system = Neutron
    group_scalar_fluxes = 'Flux_Moment_1_0_0 Flux_Moment_2_0_0 Flux_Moment_3_0_0 Flux_Moment_4_0_0'
    block = 'fuel'

    #scaling_factor = 1e0
  []
[]

# TODO: Better BC
[BCs]
  [Outer]
    type = ConvectionHeatTransferBC
    variable = T
    T_ambient = 300.0
    htc_ambient = 13.0
    boundary = '0 10000'
  []
[]

[TransportMaterials]
  [Neutronics_Fuel]
    type = FileNeutronicsMaterial
    transport_system = Neutron
    file_name = 'data/macro_xs.xml'
    source_material_id = '6'
    block = 'fuel'
    add_heating = true
  []
[]

[Materials]
  # https://www.osti.gov/servlets/purl/1330693
  [Graphite]
    type = HeatConductionMaterial
    thermal_conductivity_temperature_function = 1.0603
    temp = T
    block = 'graphite'
  []

  # https://www.osti.gov/servlets/purl/1433931
  [Fuel]
    type = HeatConductionMaterial
    thermal_conductivity_temperature_function = 0.26
    temp = T
    block = 'fuel'
  []

  # https://www.engineeringtoolbox.com/thermal-conductivity-metals-d_858.html
  [Cladding]
    type = HeatConductionMaterial
    thermal_conductivity_temperature_function = 2.36
    temp = T
    block = 'cladding'
  []

  # https://www.engineeringtoolbox.com/thermal-conductivity-metals-d_858.html
  [Stainless_Steel]
    type = HeatConductionMaterial
    thermal_conductivity_temperature_function = 0.61
    temp = T
    block = 'stainless_steel'
  []

  # https://www.engineeringtoolbox.com/air-properties-viscosity-conductivity-heat-capacity-d_1509.html
  [Air]
    type = HeatConductionMaterial
    thermal_conductivity_temperature_function = 2.624e-4
    temp = T
    block = 'air src_air'
  []
[]

[Executioner]
  type = Steady
  nl_abs_tol = 1e-11
  petsc_options_value = 'hypre boomeramg'
  petsc_options_iname = '-pc_type -pc_hypre_type'

  automatic_scaling = true
[]

[Outputs]
  exodus = true
[]
