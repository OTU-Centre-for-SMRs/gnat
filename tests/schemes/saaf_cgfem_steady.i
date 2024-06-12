[Mesh]
  [domain]
    type = CartesianMeshGenerator
    dim = 1
    dx = 10
    ix = 1000
  []
[]

[Variables]
  [angular_flux_1_1] # n = 0, Flux_Right
    order = FIRST
    family = LAGRANGE
  []
  [angular_flux_1_2] # n = 1, Flux_Left
    order = FIRST
    family = LAGRANGE
  []
[]

[AuxVariables]
  [flux_moment_1_0_0]
    order = FIRST
    family = LAGRANGE
  []
[]

[AuxKernels]
  [Flux_Moment]
    type = ParticleFluxMoment
    variable = flux_moment_1_0_0
    group_flux_ordinates = 'angular_flux_1_1 angular_flux_1_2'
    aq = AQ
    degree = 0
    order = 0
  []
[]

[Kernels]
  [Streaming_Right]
    type = SAAFStreaming
    variable = angular_flux_1_1
    aq = AQ
    ordinate_index = 0
    group_index = 0
  []
  [Streaming_Left]
    type = SAAFStreaming
    variable = angular_flux_1_2
    aq = AQ
    ordinate_index = 1
    group_index = 0
  []

  [Removal_Right]
    type = SNRemoval
    variable = angular_flux_1_1
    group_index = 0
  []
  [Removal_Left]
    type = SNRemoval
    variable = angular_flux_1_2
    group_index = 0
  []
[]

[DiracKernels]
  [Source_Right]
    type = SAAFPointSource
    variable = angular_flux_1_1
    aq = AQ
    num_groups = 1
    ordinate_index = 0
    group_index = 0
    group_source = '1000.0'
    point = '5.0 0.0 0.0'
    source_anisotropy = 0
  []
  [Source_Left]
    type = SAAFPointSource
    variable = angular_flux_1_2
    aq = AQ
    num_groups = 1
    ordinate_index = 1
    group_index = 0
    group_source = '1000.0'
    point = '5.0 0.0 0.0'
    source_anisotropy = 0
  []
[]

[BCs]
  [Right_BC]
    type = SNVacuumBC
    variable = angular_flux_1_1
    aq = AQ
    ordinate_index = 0
    boundary = 'left right'
  []
  [Left_BC]
    type = SNVacuumBC
    variable = angular_flux_1_2
    aq = AQ
    ordinate_index = 1
    boundary = 'left right'
  []
[]

[Materials]
  [Domain]
    type = AbsorbingTransportMaterial
    num_groups = 1
    group_total = 1.0
    group_speeds = 1.0
  []
[]

[UserObjects]
  [AQ]
    type = AQProvider
    dimensionality = 1D_cartesian
    aq_type = gauss_chebyshev
    n_l = 2
    n_c = 1
    major_axis = x
  []
[]

[Problem]
  type = FEProblem
[]

[Outputs]
  exodus = true
[]

[Executioner]
  type = Steady
  solve_type = PJFNK
  petsc_options_iname = '-pc_type -pc_hypre_type'
  petsc_options_value = 'hypre boomeramg'
[]
