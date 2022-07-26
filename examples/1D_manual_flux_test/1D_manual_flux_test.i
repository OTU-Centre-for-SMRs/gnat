# 1D_point_streaming.i:
# A 1D simple test case with a purely absorbing medium and a point source in the
# middle of the domain.

[Mesh]
  [domain]
    type = CartesianMeshGenerator
    dim = 1
    dx = 10
    ix = 10
  []
[]

[Variables]
  [angular_flux_1_1] # n = 0, Flux_Right
    order = FIRST
    family = MONOMIAL
  []
  [angular_flux_1_2] # n = 1, Flux_Left
    order = FIRST
    family = MONOMIAL
  []
[]

[AuxVariables]
  [flux_moment_1_0_0]
    order = FIRST
    family = MONOMIAL
  []
[]

[AuxKernels]
  [Flux_Moment]
    type = NeutronFluxMoment
    variable = flux_moment_1_0_0
    group_flux_ordinates = 'angular_flux_1_1 angular_flux_1_2'
    n_l = 2
    n_c = 1
    major_axis = x
    dimensionality = 1D_cartesian
    degree = 0
    order = 0
  []
[]

[DGKernels]
  [Upwinding_Left]
    type = ADDGNeutronStreamingUpwind
    variable = angular_flux_1_2
    n_l = 2
    n_c = 1
    major_axis = x
    dimensionality = 1D_cartesian
    ordinate_index = 1
  []
  [Upwinding_Right]
    type = ADDGNeutronStreamingUpwind
    variable = angular_flux_1_1
    n_l = 2
    n_c = 1
    major_axis = x
    dimensionality = 1D_cartesian
    ordinate_index = 0
  []
[]

[Kernels]
  [Streaming_Left]
    type = ADNeutronStreaming
    variable = angular_flux_1_2
    n_l = 2
    n_c = 1
    major_axis = x
    dimensionality = 1D_cartesian
    ordinate_index = 1
  []
  [Streaming_Right]
    type = ADNeutronStreaming
    variable = angular_flux_1_1
    n_l = 2
    n_c = 1
    major_axis = x
    dimensionality = 1D_cartesian
    ordinate_index = 0
  []

  [Removal_Left]
    type = ADNeutronRemoval
    variable = angular_flux_1_2
    group_index = 0
  []
  [Removal_Right]
    type = ADNeutronRemoval
    variable = angular_flux_1_1
    group_index = 0
  []
[]

[DiracKernels]
  [Source_Left]
    type = IsotropicNeutronPointSource
    variable = angular_flux_1_2
    n_l = 2
    n_c = 1
    major_axis = x
    dimensionality = 1D_cartesian
    intensities = '1000.0'
    points = '5.0 0.0 0.0'
  []

  [Source_Right]
    type = IsotropicNeutronPointSource
    variable = angular_flux_1_1
    n_l = 2
    n_c = 1
    major_axis = x
    dimensionality = 1D_cartesian
    intensities = '1000.0'
    points = '5.0 0.0 0.0'
  []
[]

[BCs]
  [Left_BC]
    type = ADNeutronVacuumBC
    variable = angular_flux_1_2
    n_l = 2
    n_c = 1
    major_axis = x
    dimensionality = 1D_cartesian
    ordinate_index = 1
    boundary = 'left right'
  []
  [Right_BC]
    type = ADNeutronVacuumBC
    variable = angular_flux_1_1
    n_l = 2
    n_c = 1
    major_axis = x
    dimensionality = 1D_cartesian
    ordinate_index = 0
    boundary = 'left right'
  []
[]

[Materials]
  [Domain]
    type = AbsorbingNeutronicsMaterial
    num_groups = 1
    group_removal = 1.0
    group_speeds = 1.0
  []
[]

[Problem]
  type = FEProblem
  coord_type = XYZ
[]

[Executioner]
  type = Steady
  solve_type = PJFNK
  petsc_options_iname = '-pc_type -pc_hypre_type'
  petsc_options_value = 'hypre boomeramg'
[]

[Outputs]
  exodus = true
[]
