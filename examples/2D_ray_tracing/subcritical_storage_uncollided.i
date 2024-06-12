[Mesh]
  [NeutronicsDomain]
    type = FileMeshGenerator
    file = Storage_Radiation_Mesh_Tris.e
  []
[]

[UncollidedFlux]
  [Neutron]
    uncollided_flux_treatment = ray-tracing
    num_groups = 1
    max_anisotropy = 1

    volumetric_source_blocks = fuel
    volumetric_source_moments = '40.5717'
    volumetric_source_anisotropies = '0'

    rt_n_polar = 4
  []
[]

[TransportMaterials]
  [Air]
    type = FileTransportMaterial
    transport_system = 'Neutron'
    file_name = 'macro_xs.xml'
    source_material_id = '5'
    block = air
  []
  [Cladding]
    type = FileTransportMaterial
    transport_system = 'Neutron'
    file_name = 'macro_xs.xml'
    source_material_id = '6'
    block = cladding
  []
  [Fuel]
    type = FileTransportMaterial
    transport_system = 'Neutron'
    file_name = 'macro_xs.xml'
    source_material_id = '7'
    block = fuel
  []
  [Wood]
    type = FileTransportMaterial
    transport_system = 'Neutron'
    file_name = 'macro_xs.xml'
    source_material_id = '8'
    block = wood
  []
[]

[Problem]
  type = FEProblem
  solve = false
  kernel_coverage_check = false
  material_coverage_check = false
[]

[Executioner]
  type = Steady
[]

[Outputs]
  exodus = true
  execute_on = 'TIMESTEP_END'
[]
