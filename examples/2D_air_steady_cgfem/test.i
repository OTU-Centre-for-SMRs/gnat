# A simple test case with a purely absorbing medium and a point source in the
# middle of the domain.

[Mesh]
  [gen]
    type = GeneratedMeshGenerator
    dim = 2
    xmin = 0
    xmax = 10
    ymin = 0
    ymax = 10
    nx = 101
    ny = 101
  []
[]

[TransportSystems]
  [Neutron]
    num_groups = 2
    max_anisotropy = 0
    scheme = upwinding_dfem
    particle_type = neutron
    output_angular_fluxes = true

    order = FIRST
    family = MONOMIAL

    n_azimuthal = 2
    n_polar = 2

    vacuum_boundaries = 'left right top bottom'

    point_source_locations = '5.0 5.0 0.0'
    point_source_intensities = '1000.0'
    point_source_groups = '1'
  []
[]

[TransportMaterials]
  [Domain]
    type = FileNeutronicsMaterial
    transport_system = Neutron
    file_name = './examples/2D_air_steady_cgfem/cross_sections/cross_sections.txt'
    source_material_id = '1'
  []
[]

[Problem]
  type = FEProblem
[]

[Executioner]
  type = Steady
  solve_type = NEWTON
  petsc_options_iname = '-pc_type -pc_factor_shift_type'
  petsc_options_value = 'lu       NONZERO'

  nl_abs_tol = 1e-9
  nl_max_its = 150
  line_search = 'none'

  automatic_scaling = true
  off_diagonals_in_auto_scaling = true
  compute_scaling_once = false
[]
