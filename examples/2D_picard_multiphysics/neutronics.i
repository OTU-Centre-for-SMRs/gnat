[Mesh]
  [gen]
    type = GeneratedMeshGenerator
    dim = 2
    xmin = 0
    xmax = 10
    ymin = -1
    ymax = 1
    nx = 100
    ny = 40
  []
[]

[TransportSystems]
  [Neutron]
    scheme = saaf_cfem
    particle_type = neutron
    num_groups = 2
    max_anisotropy = 0

    output_angular_fluxes = false
    disable_output = true
    flux_moment_names = 'Flux_Moment'
    scaling = 1e-14

    order = FIRST
    family = LAGRANGE

    n_azimuthal = 2
    n_polar = 2

    vacuum_boundaries = 'left right'
    source_boundaries = 'top bottom'

    boundary_source_anisotropy = '0 0'
    boundary_source_moments = '1e14 0; 1e14 0'

    constant_ic = '0.0 0.0'
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
  solve_type = NEWTON
  petsc_options_iname = '-pc_type -pc_factor_shift_type'
  petsc_options_value = 'lu       NONZERO'
  nl_abs_tol = 1e-12
  nl_max_its = 100
  line_search = 'none'

  [TimeStepper]
    type = ConstantDT
    dt = 1.0
  []
  end_time = 15
[]
