[Mesh]
  [NeutronicsDomain]
    type = FileMeshGenerator
    file = Radiation_Mesh.e
  []
[]

[TransportSystems]
  [Neutron]
    num_groups = 25
    max_anisotropy = 0
    scheme = saaf_cfem
    particle_type = neutron

    flux_moment_names = 'Flux_Moment'

    order = FIRST
    family = LAGRANGE

    n_azimuthal = 5
    n_polar = 5

    vacuum_boundaries = 'vacuum'

    volumetric_source_blocks = fuel
    volumetric_source_moments = '
    40.5717 0.0 0.0 0.0 0.0
    0.0     0.0 0.0 0.0 0.0
    0.0     0.0 0.0 0.0 0.0
    0.0     0.0 0.0 0.0 0.0
    0.0     0.0 0.0 0.0 0.0'
    volumetric_source_anisotropies = '0'
  []
[]

[TransportMaterials]
  [Air]
    type = FileNeutronicsMaterial
    transport_system = Neutron
    file_name = 'data/macro_xs.xml'
    source_material_id = '5'
    block = air
  []
  [Cladding]
    type = FileNeutronicsMaterial
    transport_system = Neutron
    file_name = 'data/macro_xs.xml'
    source_material_id = '6'
    block = cladding
  []
  [Fuel]
    type = FileNeutronicsMaterial
    transport_system = Neutron
    file_name = 'data/macro_xs.xml'
    source_material_id = '7'
    block = fuel
  []
  [Wood]
    type = FileNeutronicsMaterial
    transport_system = Neutron
    file_name = 'data/macro_xs.xml'
    source_material_id = '8'
    block = wood
  []
[]

[Executioner]
  type = Steady
  solve_type = PJFNK

  l_max_its = 50
  nl_abs_tol = 1e-12

  automatic_scaling = true
  off_diagonals_in_auto_scaling = true
  compute_scaling_once = false
[]

[Outputs]
  exodus = true
[]
