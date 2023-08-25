[Mesh]
  [NeutronicsDomain]
    type = FileMeshGenerator
    file = test_projection_1_out_Neutronics0.e
    use_for_exodus_restart = true
  []
[]

#[UserObjects]
#  [UncollidedStudy]
#    type = UncollidedFluxRayStudy
#    num_groups = 1
#
#    volumetric_source_blocks = fuel
#    volumetric_source_moments = '40.5717'
#    volumetric_source_anisotropies = '0'
#
#    n_polar = 30
#  []
#[]

#[RayKernels]
#  [UncollidedFluxKernel]
#    type = UncollidedFluxRayKernel
#    variable = UncollidedFlux
#    num_groups = 1
#    transport_system = ''
#  []
#[]

[AuxVariables]
  [UncollidedFlux]
    type = MooseVariable
    order = CONSTANT
    family = MONOMIAL
    initial_from_file_var = 'UncollidedFlux_0'
  []
[]

[Materials]
  [Air]
    type = FileNeutronicsMaterial
    transport_system = ''
    file_name = 'macro_xs.xml'
    source_material_id = '5'
    block = air
    num_groups = 1
  []
  [Cladding]
    type = FileNeutronicsMaterial
    transport_system = ''
    file_name = 'macro_xs.xml'
    source_material_id = '6'
    block = cladding
    num_groups = 1
  []
  [Fuel]
    type = FileNeutronicsMaterial
    transport_system = ''
    file_name = 'macro_xs.xml'
    source_material_id = '7'
    block = fuel
    num_groups = 1
  []
  [Wood]
    type = FileNeutronicsMaterial
    transport_system = ''
    file_name = 'macro_xs.xml'
    source_material_id = '8'
    block = wood
    num_groups = 1
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
