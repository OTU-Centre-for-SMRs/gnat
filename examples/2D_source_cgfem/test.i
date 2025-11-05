# A simple test case with a purely absorbing medium and a point source in the
# middle of the domain.

[Mesh]
  [domain]
    type = CartesianMeshGenerator
    dim = 3
    dx = '5 5'
    dy = '5 5'
    dz = '5 5'
    ix = '5 5'
    iy = '5 5'
    iz = '5 5'
    subdomain_id = '
      2 1
      1 1

      1 1
      1 1'
  []

  parallel_type = DISTRIBUTED
[]

[TransportSystems]
  [Neutron]
    scheme = saaf_cfem
    particle_type = neutron
    num_groups = 1
    max_anisotropy = 0

    order = FIRST
    family = LAGRANGE

    aq_type = level_symmetric
    ls_q_order = 18

    volumetric_source_blocks = '2'
    volumetric_source_moments = '1e0'
    volumetric_source_anisotropies = '0'
    scale_sources = true

    vacuum_boundaries = 'right top front'
    reflective_boundaries = 'left bottom back'
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
