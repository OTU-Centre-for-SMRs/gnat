[Mesh]
  [domain]
    type = CartesianMeshGenerator
    dim = 3
    dx = '10 50'
    dy = '10 90'
    dz = '10 50'
    ix = '1 5'
    iy = '1 9'
    iz = '1 5'
    subdomain_id = '
    2 1
    0 1

    1 1
    1 1'
  []
  uniform_refine = 1
[]

[TransportSystems]
  [Neutron]
    num_groups = 1
    max_anisotropy = 0
    scheme = saaf_cfem
    particle_type = neutron

    order = FIRST
    family = LAGRANGE

    n_azimuthal = 5
    n_polar = 5

    vacuum_boundaries = 'right top front'
    reflective_boundaries = 'back bottom left'

    volumetric_source_blocks = '2'
    volumetric_source_moments = '1.0'
    volumetric_source_anisotropies = '0'
  []
[]

[TransportMaterials]
  [Shield]
    type = ConstantTransportMaterial
    transport_system = 'Neutron'
    group_total = '0.1'
    group_scattering = '0.05'
    anisotropy = 0
    block = '1 2'
  []
  [Duct]
    type = ConstantTransportMaterial
    transport_system = 'Neutron'
    group_total = '1e-4'
    group_scattering = '0.5e-4'
    anisotropy = 0
    block = '0'
  []
[]

[Executioner]
  type = Steady
  solve_type = PJFNK
  petsc_options_iname = '-pc_type -pc_hypre_type -ksp_gmres_restart'
  petsc_options_value = ' hypre    boomeramg      30'
  l_max_its = 50
  nl_rel_tol = 1e-12
[]

[Outputs]
  exodus = true
  execute_on = 'TIMESTEP_END'
[]
