# The 2D Ackroyd-Riyait voided duct benchmark problem.

[Mesh]
  [domain]
    type = CartesianMeshGenerator
    dim = 2
    dx = '3 11'
    dy = '3 15'
    ix = '3 11'
    iy = '3 15'
    subdomain_id = '
    0 2
    1 2'
  []
  uniform_refine = 2
[]

[TransportSystems]
  [Neutron]
    num_groups = 1
    max_anisotropy = 0
    scheme = saaf_cfem
    particle_type = neutron

    order = FIRST
    family = LAGRANGE

    n_azimuthal = 10
    n_polar = 10

    vacuum_boundaries = 'right top'
    reflective_boundaries = 'left bottom'

    volumetric_source_blocks = '0'
    volumetric_source_moments = '1.0'
    volumetric_source_anisotropies = '0'
  []
[]

[TransportMaterials]
  [Duct]
    type = AbsorbingTransportMaterial
    transport_system = Neutron
    group_total = '0.0'
    block = '1'
  []
  [Other]
    type = AbsorbingTransportMaterial
    transport_system = Neutron
    group_total = '0.5'
    block = '0 2'
  []
[]

[Executioner]
  type = Steady
  solve_type = PJFNK
  petsc_options_iname = '-ksp_gmres_restart'
  petsc_options_value = ' 100'
  l_max_its = 50
  nl_rel_tol = 1e-12

  automatic_scaling = true
  off_diagonals_in_auto_scaling = true
  compute_scaling_once = false
[]

[Outputs]
  exodus = true
  execute_on = 'TIMESTEP_END'
[]
