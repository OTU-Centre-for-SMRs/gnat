# A simple test case with a purely absorbing medium and a point source in the
# middle of the domain.

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

[NeutronActivationStudy]
  execution_type = steady
  num_groups = 2
  max_anisotropy = 0

  [TransportSystem]
    scheme = saaf_cfem
    output_angular_fluxes = false

    order = FIRST
    family = LAGRANGE

    n_azimuthal = 2
    n_polar = 2

    vacuum_boundaries = 'left right top bottom'

    point_source_locations = '5.0 5.0 0.0'
    point_source_intensities = '1000000.0'
    point_source_groups = '1'
  []

  [NuclideSystem]
    velocity_type = constant
    constant_velocity = '0.0 4.0 0.0'

    nuclides = 'Ar40 Ar41 K41'

    order = FIRST
    family = LAGRANGE

    xs_file_name = './examples/2D_mass_transport_system_cgfem/cross_sections_micro/cross_sections.txt'
    xs_type = micro
    xs_source_material_id = '1'

    nuclide_prop_file_name = './examples/2D_mass_transport_system_cgfem/nuclide_system.txt'
    half_life_units = minutes

    density = 0.0012

    [AddNuclideBCs]
      [Ar40Inflow]
        type = ADIsotopeInflowBC
        boundary = bottom
        inflow_rate = 10.0

        excluded_nuclides = 'Ar41 K41'
      []

      [Outflow]
        type = ADIsotopeOutflowBC
        boundary = top
      []
    []
  []
[]

[Materials]
  [Domain]
    type = FileNeutronicsMaterial
    num_groups = 2
    file_name = './examples/2D_mass_transport_system_cgfem/cross_sections_macro/cross_sections.txt'
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
