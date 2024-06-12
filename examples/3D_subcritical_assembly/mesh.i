[GlobalParams]
  num_sectors_per_side = '16 16 16 16'
  background_intervals = 1
[]

[Mesh]
  [Fueled_Block]
    type = PolygonConcentricCircleMeshGenerator
    num_sides = 4
    ring_radii = '1.6383 1.8415 1.98374'
    ring_intervals = '2 1 1'
    polygon_size = 5.08

    ring_block_ids = '0 1 2'
    ring_block_names = 'fuel cladding air'
    background_block_ids = '3'
    background_block_names = 'graphite_core'

    flat_side_up = true
    quad_center_elements = true
    preserve_volumes = true

    create_outward_interface_boundaries = false
  []

  [Core_Graphite_Block]
    type = PolygonConcentricCircleMeshGenerator
    num_sides = 4
    ring_radii = '1.905 1.98374'
    ring_intervals = '2 1'
    polygon_size = 5.08

    ring_block_ids = '3 2'
    ring_block_names = 'graphite air'
    background_block_ids = '3'
    background_block_names = 'graphite_core'

    flat_side_up = true
    quad_center_elements = true
    preserve_volumes = true

    create_outward_interface_boundaries = false
  []

  [Reflector_Graphite_Block]
    type = PolygonConcentricCircleMeshGenerator
    num_sides = 4
    ring_radii = '1.905 1.98374'
    ring_intervals = '2 1'
    polygon_size = 5.08

    ring_block_ids = '4 2'
    ring_block_names = 'graphite_reflector air'
    background_block_ids = '4'
    background_block_names = 'graphite_reflector'

    flat_side_up = true
    quad_center_elements = true
    preserve_volumes = true

    create_outward_interface_boundaries = false
  []

  [Control_Block]
    type = PolygonConcentricCircleMeshGenerator
    num_sides = 4
    ring_radii = '1.905 1.98374'
    ring_intervals = '2 1'
    polygon_size = 5.08

    ring_block_ids = '5 2'
    ring_block_names = 'stainless_steel air'
    background_block_ids = '3'
    background_block_names = 'graphite_core'

    flat_side_up = true
    quad_center_elements = true
    preserve_volumes = true

    create_outward_interface_boundaries = false
  []

  [Empty_Block]
    type = PolygonConcentricCircleMeshGenerator
    num_sides = 4
    #ring_radii = '0.5 1.98374'
    #ring_intervals = '2 1'
    polygon_size = 5.08

    #ring_block_ids = '5 5'
    #ring_block_names = 'src src'
    background_block_ids = '6'
    background_block_names = 'src'

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
    #pattern_boundary = 'expanded'
    #background_block_id = '5'
    #background_block_name = 'concrete'
    #background_intervals = 3
    #square_size = 243.04

    external_boundary_id = 0
    external_boundary_name = 'vacuum'

    assign_type = 'cell'
    id_name = 'pin_id'
    generate_core_metadata = false
  []

  #[QuarterCore]
  #  type = CartesianMeshTrimmer
  #  input = 'Core_2D'
  #  center_trim_starting_index = 0
  #  center_trim_ending_index = 2
  #  center_trimming_section_boundary = 'reflective'
  #  tri_elem_subdomain_shift = 100
  #[]
  #
  [Core_3D]
    type = AdvancedExtruderGenerator
    input = 'Core_2D'
    heights = '71.12 10.16 71.12'
    num_layers = '4 2 4'
    #heights = '152.4'
    #num_layers = '10'
    direction = '0.0 0.0 1.0'

    bottom_boundary = '10001'
    top_boundary = '10000'

    subdomain_swaps = '6 2 106 102; ; 6 2 106 102'
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

#[AuxKernels]
#  [set_pin_id]
#    type = ExtraElementIDAux
#    variable = pin_id
#    extra_id_name = pin_id
#  []
#[]

[Outputs]
  exodus = true
  execute_on = timestep_end
[]
