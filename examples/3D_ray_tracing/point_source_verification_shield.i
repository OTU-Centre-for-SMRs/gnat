[Mesh]
  [domain]
    type = CartesianMeshGenerator
    dim = 3
    dx = '0.9 0.2 0.9'
    dy = '0.9 0.2 0.9'
    dz = '0.9 0.2 0.9'
    ix = '9 2 9'
    iy = '9 2 9'
    iz = '9 2 9'
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

    n_azimuthal = 1
    n_polar = 1

    vacuum_boundaries = 'left right top bottom'

    uncollided_from_multi_app = 'Uncollided'
  []
[]

[TransportMaterials]
  [Shield]
    type = ConstantTransportMaterial
    transport_system = 'Neutron'
    group_total = '1e-2'
    group_scattering = '0.0'
    anisotropy = 0
    block = '0 1'
  []
[]

[Functions]
  [Analytical]
    type = ParsedFunction
    expression = 'exp(-1e-2 * sqrt((x - 1) * (x - 1) + (y - 1) * (y - 1) + (z - 1) * (z - 1))) / (4.0 * pi * ((x - 1) * (x - 1) + (y - 1) * (y - 1) + (z - 1) * (z - 1)))'
  []
[]

[Postprocessors]
  [Error]
    type = ElementL2Error
    function = Analytical
    variable = flux_moment_1_0_0
    block = '0'
  []
  [h]
    type = AverageElementSize
    block = '0'
  []
[]

[Problem]
  solve = false
[]

[Executioner]
  type = Steady
  solve_type = PJFNK
[]

[Outputs]
  exodus = true
[]

[MultiApps]
  [Uncollided]
    type = FullSolveMultiApp
    app_type = GnatApp
    input_files = 'point_source_verification_shield_rt.i'
    execute_on = INITIAL
  []
[]
