[Mesh]
  [NeutronicsDomain]
    type = FileMeshGenerator
    file = m_and_ms.e
  []
[]

[TransportSystems]
  [Neutron]
    num_groups = 1
    max_anisotropy = 0
    scheme = saaf_cfem
    particle_type = neutron

    order = FIRST
    family = LAGRANGE

    n_azimuthal = 8
    n_polar = 8

    vacuum_boundaries = 'vacuum'

    volumetric_source_blocks = 'letters'
    volumetric_source_moments = '1.0'
    volumetric_source_anisotropies = '0'
  []
[]

[TransportMaterials]
  [Letters]
    type = ConstantTransportMaterial
    transport_system = Neutron
    group_total = '1e-1'
    group_scattering = '5e-2'
    block = 'letters'
  []
  [M_M]
    type = ConstantTransportMaterial
    transport_system = Neutron
    group_total = '5e-1'
    group_scattering = '5e-2'
    block = 'm1 m2'
  []
  [Other]
    type = ConstantTransportMaterial
    transport_system = Neutron
    group_total = '1e-4'
    group_scattering = '5e-5'
    block = 'bckgrnd'
  []
[]

[Executioner]
  type = Steady
  solve_type = PJFNK
  petsc_options_iname = '-ksp_gmres_restart'
  petsc_options_value = ' 600'
  l_max_its = 50
  nl_rel_tol = 1e-8

  automatic_scaling = true
  off_diagonals_in_auto_scaling = true
  compute_scaling_once = false
[]

[Outputs]
  exodus = true
  execute_on = 'TIMESTEP_END'
[]
