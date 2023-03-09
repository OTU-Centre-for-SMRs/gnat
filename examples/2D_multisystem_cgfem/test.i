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

[TransportSystems]
  [NeutronSN]
    scheme = saaf_cfem
    particle_type = neutron
    num_groups = 2
    max_anisotropy = 0
    output_angular_fluxes = false
    angular_flux_names = angular_flux_sn
    flux_moment_names = flux_moment_sn

    order = FIRST
    family = LAGRANGE

    n_azimuthal = 1
    n_polar = 1

    vacuum_boundaries = 'left right top bottom'

    point_source_locations = '5.0 5.0 0.0'
    point_source_intensities = '1000.0'
    point_source_groups = '1'
  []

  [NeutronDiff]
    scheme = diffusion_cfem
    particle_type = neutron
    num_groups = 2
    max_anisotropy = 0
    output_angular_fluxes = false
    angular_flux_names = angular_flux_diff
    flux_moment_names = flux_moment_diff

    order = FIRST
    family = LAGRANGE

    vacuum_boundaries = 'left right top bottom'

    point_source_locations = '5.0 5.0 0.0'
    point_source_intensities = '1000.0'
    point_source_groups = '1'
  []
[]

[TransportMaterials]
  [Water1]
    type = FileNeutronicsMaterial
    transport_system = NeutronSN
    file_name = './examples/2D_multisystem_cgfem/cross_sections/cross_sections.txt'
    source_material_id = '1'
    block = '1'
  []

  [Water2]
    type = FileNeutronicsMaterial
    transport_system = NeutronDiff
    file_name = './examples/2D_multisystem_cgfem/cross_sections/cross_sections.txt'
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
