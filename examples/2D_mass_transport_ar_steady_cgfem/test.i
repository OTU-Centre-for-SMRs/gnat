# A simple test case with a purely absorbing medium and a point source in the
# middle of the domain.

[Mesh]
  [domain]
    type = CartesianMeshGenerator
    dim = 2
    dx = 10
    dy = 10
    ix = 100
    iy = 100
  []
[]

[NeutronActivationStudy]
  num_groups = 1
  execution_type = steady
  debug_verbosity = level1

  [TransportSystem]
    scheme = saaf_cfem
    output_angular_fluxes = false

    order = FIRST
    family = LAGRANGE

    n_azimuthal = 2
    n_polar = 2

    vacuum_boundaries = 'left right top bottom'

    point_source_locations = '5.0 5.0 0.0'
    point_source_intensities = '1000.0'
    point_source_groups = '1'

    constant_ic = 0.0
  []

  [IsotopeSystem]
    velocity_type = constant
    constant_velocity = '0.0 1.0 0.0'

    isotopes = 'Ar40 Ar41'

    [AddMobileIsotopes]
      [Ar40]
        order = FIRST
        family = LAGRANGE

        isotope_name = Ar40

        diffusion_coefficient_base = 1.0
        half_life = 1.0
        half_life_units = minutes

        absorption_cross_sections = '1.0'
      []

      [Ar41]
        order = FIRST
        family = LAGRANGE

        isotope_name = Ar41
        activation_parents = Ar40

        diffusion_coefficient_base = 1.0
        half_life = 1.0
        half_life_units = minutes

        absorption_cross_sections = '0.0'
        activation_cross_sections = '1.0'
      []
    []

    [AddIsotopeBCs]
      [Inflow]
        type = ADIsotopeInflowBC
        boundary = bottom
        inflow_rate = 10.0

        excluded_isotopes = 'Ar41'
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
    type = ConstantNeutronicsMaterial
    num_groups = 1
    anisotropy = 0
    group_absorption = 1.0
    group_scattering = 1.0
    group_speeds = 220000.0
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
