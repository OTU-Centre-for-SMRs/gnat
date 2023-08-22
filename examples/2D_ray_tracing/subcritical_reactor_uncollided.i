# TODO: Fix the mesh. QUADSHELL4 is not supported in raytracing.
# Maybe add QUADSHELL4 as a raytracable element?
[Mesh]
  [NeutronicsDomain]
    type = FileMeshGenerator
    file = Reactor_Radiation_Mesh.e
  []
[]

[UserObjects]
  [UncollidedStudy]
    type = UncollidedFluxRayStudy
    num_groups = 1
    group_index = 0

    volumetric_source_blocks = 3
    volumetric_source_moments = '40.5717'
    volumetric_source_anisotropies = '0'

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
  [Air]
    type = FileNeutronicsMaterial
    transport_system = ''
    num_groups = 1
    file_name = 'xs_macro_reactor/air_cross_sections.txt'
    source_material_id = '5' #
    block = '2 7 8 9 '
  []

  [Cladding]
    type = FileNeutronicsMaterial
    transport_system = ''
    num_groups = 1
    file_name = 'xs_macro_reactor/clad_cross_sections.txt'
    source_material_id = '7' #
    block = 4
  []

  [Fuel]
    type = FileNeutronicsMaterial
    transport_system = ''
    num_groups = 1
    file_name = 'xs_macro_reactor/fuel_cross_sections.txt'
    source_material_id = '6' #
    block = 3
  []

  [Graphite]
    type = FileNeutronicsMaterial
    transport_system = ''
    num_groups = 1
    file_name = 'xs_macro_reactor/mod_cross_sections.txt'
    source_material_id = '8' #
    block = '10 11'
  []

  [Control]
    type = FileNeutronicsMaterial
    transport_system = ''
    num_groups = 1
    file_name = 'xs_macro_reactor/steel_cross_sections.txt'
    source_material_id = '9' #
    block = 6
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
