# A simple test case with a purely absorbing medium and a point source in the
# middle of the domain.

[Mesh]
  [domain]
    type = CartesianMeshGenerator
    dim = 3
    dx = '10 20 10 20'
    dy = '10 40 10 40'
    dz = '10 20 10 20'
    ix = '1 2 1 2'
    iy = '1 4 1 4'
    iz = '1 2 1 2'
    subdomain_id = '
    2 1 1 1
    0 1 1 1
    0 0 0 1
    1 1 1 1

    1 1 1 1
    1 1 1 1
    1 1 0 1
    1 1 1 1

    1 1 1 1
    1 1 1 1
    1 1 0 1
    1 1 0 1

    1 1 1 1
    1 1 1 1
    1 1 1 1
    1 1 1 1'
  []
  uniform_refine = 4
  #parallel_type = DISTRIBUTED
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

[Materials]
  [Duct]
    type = AbsorbingNeutronicsMaterial
    transport_system = 'Neutron'
    num_groups = 1
    group_total = '1e-4'
    block = 0
  []
  [Shield]
    type = AbsorbingNeutronicsMaterial
    transport_system = 'Neutron'
    num_groups = 1
    group_total = '0.1'
    block = '1 2'
  []
[]

[Problem]
  solve = false
  kernel_coverage_check = false
[]

[Executioner]
  type = Steady
[]

[Outputs]
  exodus = true
  execute_on = TIMESTEP_END
[]
