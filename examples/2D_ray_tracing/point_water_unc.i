# A simple test case with a purely absorbing medium and a point source in the
# middle of the domain.

[Mesh]
  [domain]
    type = CartesianMeshGenerator
    dim = 2
    dx = '5.0 10.0'
    dy = '10.0 5.0'
    ix = '5 10'
    iy = '10 5'
    subdomain_id = '
      1 1
      2 1'
  []
  uniform_refine = 2
[]

[UncollidedFlux]
  [Neutron]
    uncollided_flux_treatment = ray-tracing
    num_groups = 2
    max_anisotropy = 0

    volumetric_source_blocks = '2'
    volumetric_source_moments = '1.0 0.5'
    volumetric_source_anisotropies = '0'
  []
[]

[TransportMaterials]
  [Water]
    type = FileNeutronicsMaterial
    transport_system = Neutron
    file_name = 'macro_xs_water.xml'
    source_material_id = '1'
    block = '1'
  []
  [Air]
    type = FileNeutronicsMaterial
    transport_system = Neutron
    file_name = 'macro_xs_air.xml'
    source_material_id = '5'
    block = '2'
  []
[]

[Problem]
  type = FEProblem
  kernel_coverage_check = false
  solve = false
[]

[Executioner]
  type = Steady
  solve_type = PJFNK
[]

[Outputs]
  exodus = true
[]
