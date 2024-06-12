# A simple test case with a purely absorbing medium and a point source in the
# middle of the domain.

[Mesh]
  [domain]
    type = CartesianMeshGenerator
    dim = 2
    dx = 10
    dy = 10
    ix = 101
    iy = 101
  []
[]

[TransportSystems]
  [Sub_Neutron]
    num_groups = 1
    scheme = saaf_cfem
    particle_type = neutron
    flux_moment_names = sub_flux_moment

    order = FIRST
    family = LAGRANGE

    n_azimuthal = 2
    n_polar = 2

    max_anisotropy = 0
    vacuum_boundaries = 'left right top bottom'

    point_source_locations = '5.0 5.0 0.0'
    point_source_moments = '1000.0'
    point_source_anisotropies = '0'
    scale_sources = true
  []
[]

[TransportMaterials]
  [Domain]
    type = AbsorbingTransportMaterial
    transport_system = Sub_Neutron
    group_total = 2.0
    group_speeds = 2200.0
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
