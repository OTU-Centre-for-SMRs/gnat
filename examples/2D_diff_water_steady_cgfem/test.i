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
      2 2 2
      2 1 1
      2 2 2'
  []
[]

[TransportSystems]
  [Neutron]
    num_groups = 2
    scheme = diffusion_cfem
    particle_type = neutron

    order = FIRST
    family = LAGRANGE

    vacuum_boundaries = 'left right top bottom'

    point_source_locations = '5.0 5.0 0.0'
    point_source_moments = '1e3 0.0'
    point_source_anisotropies = '0'
    scale_sources = true
  []
[]

[TransportMaterials]
  [Water]
    type = FileNeutronicsMaterial
    transport_system = Neutron
    file_name = 'cross_sections/cross_sections.txt'
    source_material_id = '1'
    block = '1'
  []

  [Shield]
    type = AbsorbingNeutronicsMaterial
    transport_system = Neutron
    group_absorption = '1.0 3.0'
    group_speeds = '10000000000.0 2200000.0'
    block = '2'
  []
[]

[Executioner]
  type = Steady
  solve_type = PJFNK
  petsc_options_iname = '-pc_type -pc_hypre_type -ksp_gmres_restart'
  petsc_options_value = 'hypre boomeramg 10'

  l_max_its = 50
  nl_rel_tol = 1e-12
[]

[Outputs]
  exodus = true
[]
