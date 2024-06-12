# A simple test case with a purely absorbing medium and a point source in the
# middle of the domain.

[Mesh]
  [domain]
    type = CartesianMeshGenerator
    dim = 2
    dx = '4 20 20'
    dy = '4 20 20'
    ix = '4 20 20'
    iy = '4 20 20'
    subdomain_id = '
      2 1 1
      1 1 1
      1 1 1'
  []
[]

[TransportSystems]
  [Neutron]
    scheme = saaf_cfem
    particle_type = neutron
    num_groups = 1
    max_anisotropy = 0

    order = FIRST
    family = LAGRANGE

    n_azimuthal = 20
    n_polar = 20

    volumetric_source_blocks = '2'
    volumetric_source_moments = '1e0'
    volumetric_source_anisotropies = '0'
    scale_sources = true

    vacuum_boundaries = 'left right top bottom'
  []
[]

[TransportMaterials]
  [Domain]
    type = ConstantTransportMaterial
    transport_system = Neutron
    group_total = '1.1e-6'
    group_scattering = '1e-6'
    block = '1 2'
  []
[]

[Executioner]
  type = Steady
  solve_type = PJFNK
  #petsc_options_iname = '-pc_type -pc_hypre_type -ksp_gmres_restart'
  #petsc_options_value = ' hypre    boomeramg      600'
  petsc_options_iname = '-ksp_gmres_restart'
  petsc_options_value = ' 600'
  l_max_its = 50
  nl_rel_tol = 1e-12
[]

[Outputs]
  exodus = true
  execute_on = 'TIMESTEP_END'
[]
