[Mesh]
  [domain]
    type = CartesianMeshGenerator
    dim = 1
    dx = 10
    ix = 100
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
  [Upwinding_Right]
    type = ADDFEMUpwinding
    variable = angular_flux_1_1
    n_l = 2
    n_c = 1
    major_axis = x
    dimensionality = 1D_cartesian
    ordinate_index = 0
  []
  [Upwinding_Left]
    type = ADDFEMUpwinding
    variable = angular_flux_1_2
    n_l = 2
    n_c = 1
    major_axis = x
    dimensionality = 1D_cartesian
    ordinate_index = 1
  []
[]

[Kernels]
  [Streaming_Right]
    type = ADDFEMStreaming
    variable = angular_flux_1_1
    n_l = 2
    n_c = 1
    major_axis = x
    dimensionality = 1D_cartesian
    ordinate_index = 0
  []
  [Streaming_Left]
    type = ADDFEMStreaming
    variable = angular_flux_1_2
    n_l = 2
    n_c = 1
    major_axis = x
    dimensionality = 1D_cartesian
    ordinate_index = 1
  []

  [Removal_Right]
    type = ADSNRemoval
    variable = angular_flux_1_1
    group_index = 0
  []
  [Removal_Left]
    type = ADSNRemoval
    variable = angular_flux_1_2
    group_index = 0
  []
[]

[DiracKernels]
  [Source_Right]
    type = DFEMIsoPointSource
    variable = angular_flux_1_1
    n_l = 2
    n_c = 1
    major_axis = x
    dimensionality = 1D_cartesian
    intensities = '1000.0'
    points = '5.0 0.0 0.0'
  []
  [Source_Left]
    type = DFEMIsoPointSource
    variable = angular_flux_1_2
    n_l = 2
    n_c = 1
    major_axis = x
    dimensionality = 1D_cartesian
    intensities = '1000.0'
    points = '5.0 0.0 0.0'
  []
[]

[BCs]
  [Right_BC]
    type = ADSNVacuumBC
    variable = angular_flux_1_1
    n_l = 2
    n_c = 1
    major_axis = x
    dimensionality = 1D_cartesian
    ordinate_index = 0
    boundary = 'left right'
  []
  [Left_BC]
    type = ADSNVacuumBC
    variable = angular_flux_1_2
    n_l = 2
    n_c = 1
    major_axis = x
    dimensionality = 1D_cartesian
    ordinate_index = 1
    boundary = 'left right'
  []
[]

[Materials]
  [Domain]
    type = AbsorbingNeutronicsMaterial
    num_groups = 1
    group_absorption = 0.0
    group_speeds = 1.0
  []
[]

[Problem]
  type = FEProblem
  coord_type = XYZ
[]

[Executioner]
  type = Steady
  solve_type = NEWTON
[]

[Outputs]
  exodus = true
[]
