[Mesh]
  [Result]
    type = FileMeshGenerator
    file = 'neutronics_rod_g26_L0_corrected_out.e'
    use_for_exodus_restart = true
  []
[]

[AuxVariables]
  [G1]
    type = MooseVariable
    order = FIRST
    family = LAGRANGE
    initial_from_file_var = 'Flux_Moment_1_0_0'
  []

  [G2]
    type = MooseVariable
    order = FIRST
    family = LAGRANGE
    initial_from_file_var = 'Flux_Moment_2_0_0'
  []

  [G3]
    type = MooseVariable
    order = FIRST
    family = LAGRANGE
    initial_from_file_var = 'Flux_Moment_3_0_0'
  []

  [G4]
    type = MooseVariable
    order = FIRST
    family = LAGRANGE
    initial_from_file_var = 'Flux_Moment_4_0_0'
  []

  [G5]
    type = MooseVariable
    order = FIRST
    family = LAGRANGE
    initial_from_file_var = 'Flux_Moment_5_0_0'
  []

  [G6]
    type = MooseVariable
    order = FIRST
    family = LAGRANGE
    initial_from_file_var = 'Flux_Moment_6_0_0'
  []

  [G7]
    type = MooseVariable
    order = FIRST
    family = LAGRANGE
    initial_from_file_var = 'Flux_Moment_7_0_0'
  []

  [G8]
    type = MooseVariable
    order = FIRST
    family = LAGRANGE
    initial_from_file_var = 'Flux_Moment_8_0_0'
  []

  [G9]
    type = MooseVariable
    order = FIRST
    family = LAGRANGE
    initial_from_file_var = 'Flux_Moment_9_0_0'
  []

  [G10]
    type = MooseVariable
    order = FIRST
    family = LAGRANGE
    initial_from_file_var = 'Flux_Moment_10_0_0'
  []

  [G11]
    type = MooseVariable
    order = FIRST
    family = LAGRANGE
    initial_from_file_var = 'Flux_Moment_11_0_0'
  []

  [G12]
    type = MooseVariable
    order = FIRST
    family = LAGRANGE
    initial_from_file_var = 'Flux_Moment_12_0_0'
  []
  [G13]
    type = MooseVariable
    order = FIRST
    family = LAGRANGE
    initial_from_file_var = 'Flux_Moment_13_0_0'
  []

  [G14]
    type = MooseVariable
    order = FIRST
    family = LAGRANGE
    initial_from_file_var = 'Flux_Moment_14_0_0'
  []

  [G15]
    type = MooseVariable
    order = FIRST
    family = LAGRANGE
    initial_from_file_var = 'Flux_Moment_15_0_0'
  []

  [G16]
    type = MooseVariable
    order = FIRST
    family = LAGRANGE
    initial_from_file_var = 'Flux_Moment_16_0_0'
  []

  [G17]
    type = MooseVariable
    order = FIRST
    family = LAGRANGE
    initial_from_file_var = 'Flux_Moment_17_0_0'
  []

  [G18]
    type = MooseVariable
    order = FIRST
    family = LAGRANGE
    initial_from_file_var = 'Flux_Moment_18_0_0'
  []

  [G19]
    type = MooseVariable
    order = FIRST
    family = LAGRANGE
    initial_from_file_var = 'Flux_Moment_19_0_0'
  []

  [G20]
    type = MooseVariable
    order = FIRST
    family = LAGRANGE
    initial_from_file_var = 'Flux_Moment_20_0_0'
  []

  [G21]
    type = MooseVariable
    order = FIRST
    family = LAGRANGE
    initial_from_file_var = 'Flux_Moment_21_0_0'
  []

  [G22]
    type = MooseVariable
    order = FIRST
    family = LAGRANGE
    initial_from_file_var = 'Flux_Moment_22_0_0'
  []

  [G23]
    type = MooseVariable
    order = FIRST
    family = LAGRANGE
    initial_from_file_var = 'Flux_Moment_23_0_0'
  []

  [G24]
    type = MooseVariable
    order = FIRST
    family = LAGRANGE
    initial_from_file_var = 'Flux_Moment_24_0_0'
  []

  [G25]
    type = MooseVariable
    order = FIRST
    family = LAGRANGE
    initial_from_file_var = 'Flux_Moment_25_0_0'
  []

  [G26]
    type = MooseVariable
    order = FIRST
    family = LAGRANGE
    initial_from_file_var = 'Flux_Moment_26_0_0'
  []
[]

[Postprocessors]
  [G1_Int]
    type = ElementIntegralVariablePostprocessor
    variable = G1
    #block = '5'
  []

  [G2_Int]
    type = ElementIntegralVariablePostprocessor
    variable = G2
    #block = '5'
  []

  [G3_Int]
    type = ElementIntegralVariablePostprocessor
    variable = G3
    #block = '5'
  []

  [G4_Int]
    type = ElementIntegralVariablePostprocessor
    variable = G4
    #block = '5'
  []

  [G5_Int]
    type = ElementIntegralVariablePostprocessor
    variable = G5
    #block = '5'
  []

  [G6_Int]
    type = ElementIntegralVariablePostprocessor
    variable = G6
    #block = '5'
  []

  [G7_Int]
    type = ElementIntegralVariablePostprocessor
    variable = G7
    #block = '5'
  []

  [G8_Int]
    type = ElementIntegralVariablePostprocessor
    variable = G8
    #block = '5'
  []

  [G9_Int]
    type = ElementIntegralVariablePostprocessor
    variable = G9
    #block = '5'
  []

  [G10_Int]
    type = ElementIntegralVariablePostprocessor
    variable = G10
    #block = '5'
  []

  [G11_Int]
    type = ElementIntegralVariablePostprocessor
    variable = G11
    #block = '5'
  []

  [G12_Int]
    type = ElementIntegralVariablePostprocessor
    variable = G12
    #block = '5'
  []

  [G13_Int]
    type = ElementIntegralVariablePostprocessor
    variable = G13
    #block = '5'
  []

  [G14_Int]
    type = ElementIntegralVariablePostprocessor
    variable = G14
    #block = '5'
  []

  [G15_Int]
    type = ElementIntegralVariablePostprocessor
    variable = G15
    #block = '5'
  []

  [G16_Int]
    type = ElementIntegralVariablePostprocessor
    variable = G16
    #block = '5'
  []

  [G17_Int]
    type = ElementIntegralVariablePostprocessor
    variable = G17
    #block = '5'
  []

  [G18_Int]
    type = ElementIntegralVariablePostprocessor
    variable = G18
    #block = '5'
  []

  [G19_Int]
    type = ElementIntegralVariablePostprocessor
    variable = G19
    #block = '5'
  []

  [G20_Int]
    type = ElementIntegralVariablePostprocessor
    variable = G20
    #block = '5'
  []

  [G21_Int]
    type = ElementIntegralVariablePostprocessor
    variable = G21
    #block = '5'
  []

  [G22_Int]
    type = ElementIntegralVariablePostprocessor
    variable = G22
    #block = '5'
  []

  [G23_Int]
    type = ElementIntegralVariablePostprocessor
    variable = G23
    #block = '5'
  []

  [G24_Int]
    type = ElementIntegralVariablePostprocessor
    variable = G24
    #block = '5'
  []

  [G25_Int]
    type = ElementIntegralVariablePostprocessor
    variable = G25
    #block = '5'
  []

  [G26_Int]
    type = ElementIntegralVariablePostprocessor
    variable = G26
    #block = '5'
  []
[]

[Outputs]
  csv = true
[]

[Executioner]
  type = Steady
[]

[Problem]
  solve = false
  kernel_coverage_check = false
  material_coverage_check = false
[]

