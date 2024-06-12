# The 2D Ackroyd-Riyait voided dog-legged duct benchmark problem.

[Mesh]
  [NeutronicsDomain]
    type = FileMeshGenerator
    file = kirk_name.e
  []
  uniform_refine = 0
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

    vacuum_boundaries = 'vacuum'

    volumetric_source_blocks = 'src'
    volumetric_source_moments = '1.0'
    volumetric_source_anisotropies = '0'
  []
[]

[TransportMaterials]
  [Kirk]
    type = ConstantTransportMaterial
    transport_system = Neutron
    group_total = '1e-1'
    group_scattering = '5e-2'
    block = 'src'
  []
  [Other]
    type = ConstantTransportMaterial
    transport_system = Neutron
    group_total = '1e-1'
    group_scattering = '5e-2'
    block = 'shield'
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
