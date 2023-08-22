[Mesh]
  [NeutronicsDomain]
    type = FileMeshGenerator
    file = Storage_Radiation_Mesh_Tris.e
  []
  uniform_refine = 0
[]

[UserObjects]
  [UncollidedStudy]
    type = UncollidedFluxRayStudy
    num_groups = 1
    group_index = 0

    volumetric_source_blocks = fuel
    volumetric_source_moments = '40.5717'
    volumetric_source_anisotropies = '0'

    n_polar = 5
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
  [Air]
    type = FileNeutronicsMaterial
    transport_system = ''
    file_name = 'xs_macro_storage/air_cross_sections.txt'
    source_material_id = '5'
    block = air
    num_groups = 1
  []
  [Cladding]
    type = FileNeutronicsMaterial
    transport_system = ''
    file_name = 'xs_macro_storage/clad_cross_sections.txt'
    source_material_id = '6'
    block = cladding
    num_groups = 1
  []
  [Fuel]
    type = FileNeutronicsMaterial
    transport_system = ''
    file_name = 'xs_macro_storage/fuel_cross_sections.txt'
    source_material_id = '7'
    block = fuel
    num_groups = 1
  []
  [Wood]
    type = FileNeutronicsMaterial
    transport_system = ''
    file_name = 'xs_macro_storage/box_cross_sections.txt'
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
