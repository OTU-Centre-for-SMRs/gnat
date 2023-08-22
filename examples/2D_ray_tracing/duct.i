[Mesh]
  [domain]
    type = CartesianMeshGenerator
    dim = 2
    dx = '11 3 3 11'
    dy = '15 3 3 15'
    ix = '11 3 3 11'
    iy = '15 3 3 15'
    subdomain_id = '
    1 0 0 1
    1 2 2 1
    1 2 2 1
    1 0 0 1'
  []
  uniform_refine = 1
[]

[UserObjects]
  [UncollidedStudy]
    type = UncollidedFluxRayStudy
    num_groups = 1
    group_index = 0

    point_source_locations = '
    14.0  0.0 0.0
    14.0 36.0 0.0'
    point_source_moments = '
    1.0;
    1.0'
    point_source_anisotropies = '0 0'

    source_boundaries = '1 3'
    boundary_source_anisotropy = '0 0'
    boundary_source_moments = '1.0; 1.0'

    volumetric_source_blocks = '2'
    volumetric_source_anisotropies = '0'
    volumetric_source_moments = '1.0'

    n_polar = 30
  []
[]

[RayKernels]
  [UncollidedFluxKernel]
    type = UncollidedFluxRayKernel
    variable = UncollidedFlux
    group_index = 0
    transport_system = ''
  []
[]

[AuxVariables]
  [UncollidedFlux]
    type = MooseVariable
    order = CONSTANT
    family = MONOMIAL
  []
[]

[Materials]
  [Shield]
    type = AbsorbingNeutronicsMaterial
    transport_system = ''
    num_groups = 1
    group_absorption = '0.5'
    group_speeds = '220000'
    block = '1 2'
  []
  [Duct]
    type = AbsorbingNeutronicsMaterial
    transport_system = ''
    num_groups = 1
    group_absorption = '0.0'
    group_speeds = '220000'
    block = '0'
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
[]
