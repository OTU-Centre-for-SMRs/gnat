[Mesh]
  [domain]
    type = FileMeshGenerator
    file = 'reference_in.e'
  []
  #[domain]
  #  type = CartesianMeshGenerator
  #  dim = 2
  #  dx = '1.0 5.0 1.0 3.0'
  #  dy = '1.0 5.0 1.0 3.0'
  #  ix = '2 10 2 6'
  #  iy = '2 10 2 6'
  #  subdomain_id = '2 0 0 0
  #                  0 0 0 0
  #                  0 0 1 0
  #                  0 0 0 0'
  #[]
  #uniform_refine = 5
[]

[UncollidedFlux]
  [Neutron]
    uncollided_flux_treatment = ray-tracing
    num_groups = 1
    max_anisotropy = 0

    point_source_locations = '0.0 0.0 0.0'
    point_source_moments = '1.0'
    point_source_anisotropies = '0'
  []
[]

[TransportMaterials]
  [Shield]
    type = AbsorbingTransportMaterial
    transport_system = 'Neutron'
    #group_total = '1.0'
    group_total = '10.0'
    block = '1'
  []
  [Empty]
    type = AbsorbingTransportMaterial
    transport_system = 'Neutron'
    group_total = '0.0'
    block = '0 2'
  []
[]

[Problem]
  type = FEProblem
  solve = false
  kernel_coverage_check = false
[]

[Executioner]
  type = Steady
[]

[Outputs]
  exodus = true
  execute_on = 'TIMESTEP_END'
[]
