# Coupled mass transport with neutron activation.

[Mesh]
  [domain]
    type = CartesianMeshGenerator
    dim = 2
    dx = '3.5 3 3.5'
    dy = '3.5 3 3.5'
    ix = '35 31 35'
    iy = '35 31 35'
    subdomain_id = '
      1 1 1
      1 1 1
      1 1 1'
  []
[]

[TransportSystems]
  [Neutron]
    scheme = saaf_cfem
    particle_type = neutron
    num_groups = 2
    max_anisotropy = 0

    order = FIRST
    family = LAGRANGE

    n_azimuthal = 2
    n_polar = 2

    vacuum_boundaries = 'left right top bottom'

    point_source_locations = '5.0 5.0 0.0'
    point_source_intensities = '1000.0'
    point_source_groups = '1'
  []
[]

[NuclideSystem]
  velocity_type = constant
  constant_velocity = '0.0 1.0 0.0'

  transport_system = Neutron

  nuclides = 'Ar40 Ar41 K41'

  order = FIRST
  family = LAGRANGE

  xs_file_name = './examples/2D_mass_transport_system_steady_cgfem/xs_micro/cross_sections.txt'
  xs_type = micro
  xs_source_material_id = '1'

  nuclide_prop_file_name = './examples/2D_mass_transport_system_steady_cgfem/nuclide_system.txt'
  half_life_units = minutes

  density = 0.0012

  [AddNuclideBCs]
    [Inflow]
      type = ADIsotopeInflowBC
      boundary = bottom
      inflow_rate = 0.000011163
      # atom fraction Ar40 * weight fraction Ar * density of air = Inflow mass concentration of Ar40
      # 0.996035 * 0.009340 * 0.0012 = 0.000011163

      excluded_nuclides = 'Ar41 K41'
    []

    [Outflow]
      type = ADIsotopeOutflowBC
      boundary = top
    []
  []
[]

[TransportMaterials]
  [Domain]
    type = FileNeutronicsMaterial
    transport_system = Neutron
    file_name = './examples/2D_mass_transport_system_steady_cgfem/xs_macro/cross_sections.txt'
    source_material_id = '1'
    block = '1'
  []
[]

[Problem]
  type = FEProblem
[]

[Executioner]
  type = Steady
  solve_type = PJFNK
  petsc_options_iname = '-pc_type -pc_hypre_type -ksp_gmres_restart'
  petsc_options_value = 'hypre boomeramg 10'
  l_max_its = 50
  nl_rel_tol = 1e-12
[]
