[Mesh]
  [domain]
    type = CartesianMeshGenerator
    dim = 3
    dx = '0.5 1.0 0.5'
    dy = '0.5 1.0 0.5'
    dz = '0.5 1.0 0.5'
    ix = '5 10 5'
    iy = '5 10 5'
    iz = '5 10 5'
    subdomain_id = '
    0 0 0
    0 0 0
    0 0 0

    0 0 0
    0 1 0
    0 0 0

    0 0 0
    0 0 0
    0 0 0'
  []
  uniform_refine = 3
[]

[UncollidedFlux]
  [Neutron]
    uncollided_flux_treatment = ray-tracing
    num_groups = 1
    max_anisotropy = 0

    point_source_locations = '1.0 1.0 1.0'
    point_source_moments = '1.0'
    point_source_anisotropies = '0'
  []
[]

[TransportMaterials]
  [Duct]
    type = ConstantTransportMaterial
    transport_system = 'Neutron'
    group_total = '0.0'
    group_scattering = '0.0'
    anisotropy = 0
    block = '0 1'
  []
[]

[Functions]
  [Analytical]
    type = ParsedFunction
    expression = '1.0 / (4.0 * pi * ((x - 1) * (x - 1) + (y - 1) * (y - 1) + (z - 1) * (z - 1)))'
  []
[]

[Postprocessors]
  [Error]
    type = ElementL2Error
    function = Analytical
    variable = uncollided_flux_moment_1_0_0
    block = '0'
  []
  [h]
    type = AverageElementSize
    block = '0'
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
[]
