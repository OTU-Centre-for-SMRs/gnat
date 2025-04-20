#----------------------------------------------------------------------------------------
# Core geometrical information
#----------------------------------------------------------------------------------------
fuel_radius          = 1.6383
clad_thickness       = 0.2032
channel_radius       = 1.98374
graphite_plug_radius = 1.905
control_rod_radius   = 1.905
graphite_block_lw    = 10.16
core_height          = ${fparse 2.0 * 76.2}
#----------------------------------------------------------------------------------------

#----------------------------------------------------------------------------------------
# Meshing parameters
#----------------------------------------------------------------------------------------
NUM_SECTORS              = 4
FUEL_RADIAL_DIVISIONS    = 4
CONTROL_RADIAL_DIVISIONS = 4
BACKGROUND_DIVISIONS     = 2
AXIAL_DIVISIONS          = 5
#----------------------------------------------------------------------------------------

[Mesh]
  [Fueled_Block]
    type = PolygonConcentricCircleMeshGenerator
    num_sides = 4
    num_sectors_per_side = '${NUM_SECTORS} ${NUM_SECTORS} ${NUM_SECTORS} ${NUM_SECTORS}'
    ring_radii = '${fuel_radius} ${fparse fuel_radius + clad_thickness} ${channel_radius}'
    ring_intervals = '${FUEL_RADIAL_DIVISIONS} 1 1'
    polygon_size = ${fparse graphite_block_lw / 2.0}

    ring_block_ids = '0 1 2'
    ring_block_names = 'fuel cladding air'
    background_block_ids = '3'
    background_block_names = 'graphite_core'
    background_intervals = ${BACKGROUND_DIVISIONS}

    flat_side_up = true
    quad_center_elements = true
    preserve_volumes = true

    create_outward_interface_boundaries = false
  []

  [Core_Graphite_Block]
    type = PolygonConcentricCircleMeshGenerator
    num_sides = 4
    num_sectors_per_side = '${NUM_SECTORS} ${NUM_SECTORS} ${NUM_SECTORS} ${NUM_SECTORS}'
    ring_radii = '${graphite_plug_radius} ${channel_radius}'
    ring_intervals = '2 1'
    polygon_size = ${fparse graphite_block_lw / 2.0}

    ring_block_ids = '3 2'
    ring_block_names = 'graphite air'
    background_block_ids = '3'
    background_block_names = 'graphite_core'
    background_intervals = ${BACKGROUND_DIVISIONS}

    flat_side_up = true
    quad_center_elements = true
    preserve_volumes = true

    create_outward_interface_boundaries = false
  []

  [Reflector_Graphite_Block]
    type = PolygonConcentricCircleMeshGenerator
    num_sides = 4
    num_sectors_per_side = '${NUM_SECTORS} ${NUM_SECTORS} ${NUM_SECTORS} ${NUM_SECTORS}'
    ring_radii = '${graphite_plug_radius} ${channel_radius}'
    ring_intervals = '2 1'
    polygon_size = ${fparse graphite_block_lw / 2.0}

    ring_block_ids = '4 2'
    ring_block_names = 'graphite_reflector air'
    background_block_ids = '4'
    background_block_names = 'graphite_reflector'
    background_intervals = ${BACKGROUND_DIVISIONS}

    flat_side_up = true
    quad_center_elements = true
    preserve_volumes = true

    create_outward_interface_boundaries = false
  []

  [Control_Block]
    type = PolygonConcentricCircleMeshGenerator
    num_sides = 4
    num_sectors_per_side = '${NUM_SECTORS} ${NUM_SECTORS} ${NUM_SECTORS} ${NUM_SECTORS}'
    ring_radii = '${control_rod_radius} ${channel_radius}'
    ring_intervals = '${CONTROL_RADIAL_DIVISIONS} 1'
    polygon_size = ${fparse graphite_block_lw / 2.0}

    ring_block_ids = '5 2'
    ring_block_names = 'stainless_steel air'
    background_block_ids = '3'
    background_block_names = 'graphite_core'
    background_intervals = ${BACKGROUND_DIVISIONS}

    flat_side_up = true
    quad_center_elements = true
    preserve_volumes = true

    create_outward_interface_boundaries = false
  []

  [Empty_Block]
    type = PolygonConcentricCircleMeshGenerator
    num_sides = 4
    num_sectors_per_side = '${NUM_SECTORS} ${NUM_SECTORS} ${NUM_SECTORS} ${NUM_SECTORS}'
    polygon_size = ${fparse graphite_block_lw / 2.0}

    background_block_ids = '6'
    background_block_names = 'src'
    background_intervals = ${BACKGROUND_DIVISIONS}

    flat_side_up = true
    quad_center_elements = true
    preserve_volumes = true

    create_outward_interface_boundaries = false
  []

  [Core_2D]
    type = PatternedCartesianMeshGenerator
    inputs = 'Fueled_Block Reflector_Graphite_Block Core_Graphite_Block Control_Block Empty_Block'
    pattern = ' 1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1;
                1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1;
                1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1;
                1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1;
                1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1;
                1  1  1  1  1  0  2  0  2  0  2  0  2  0  1  1  1  1  1;
                1  1  1  1  1  2  0  2  0  2  0  2  0  2  1  1  1  1  1;
                1  1  1  1  1  0  2  0  2  0  2  0  2  0  1  1  1  1  1;
                1  1  1  1  1  2  0  2  2  3  2  2  0  2  1  1  1  1  1;
                1  1  1  1  1  0  2  0  3  4  3  0  2  0  1  1  1  1  1;
                1  1  1  1  1  2  0  2  2  3  2  2  0  2  1  1  1  1  1;
                1  1  1  1  1  0  2  0  2  0  2  0  2  0  1  1  1  1  1;
                1  1  1  1  1  2  0  2  0  2  0  2  0  2  1  1  1  1  1;
                1  1  1  1  1  0  2  0  2  0  2  0  2  0  1  1  1  1  1;
                1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1;
                1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1;
                1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1;
                1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1;
                1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1'

    pattern_boundary = 'none'

    external_boundary_id = 0
    external_boundary_name = 'vacuum'

    assign_type = 'cell'
    id_name = 'pin_id'
    generate_core_metadata = false
  []

  [2D_Trimmed]
    type = CartesianMeshTrimmer
    input = Core_2D
    center_trim_starting_index = 0
    center_trim_ending_index = 2
    tri_elem_subdomain_shift = '1000'
    center_trimming_section_boundary = 'reflective'
  []

  #[Delete_Air]
  #  type = BlockDeletionGenerator
  #  input = '2D_Trimmed'
  #  block = 'air src src_trimmer_tri'
  #[]

  [Core_3D]
    type = AdvancedExtruderGenerator
    input = '2D_Trimmed'
    heights = '${fparse core_height / 2.0}'
    num_layers = '${AXIAL_DIVISIONS}'
    direction = '0.0 0.0 1.0'

    bottom_boundary = '10001'
    top_boundary = '10000'
  []
[]

[Executioner]
  type = Steady
[]

[Problem]
  solve = false
[]

[AuxVariables]
  [pin_id]
    family = MONOMIAL
    order = CONSTANT
  []
[]

[AuxKernels]
  [set_pin_id]
    type = ExtraElementIDAux
    variable = pin_id
    extra_id_name = pin_id
  []
[]

[Outputs]
  exodus = true
  execute_on = timestep_end
[]
