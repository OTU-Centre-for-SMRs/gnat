[Mesh]
  [bottom]
    type = CartesianMeshGenerator
    dim = 2
    dx = '9'
    dy = '1'
    ix = '18'
    iy = '2'
  []
  [bottom_trans]
    type = TransformGenerator
    input = bottom
    transform = TRANSLATE
    vector_value = '1.0 0.0 0.0'
  []
  [bottom_rename]
    type = RenameBoundaryGenerator
    input = bottom_trans
    old_boundary = 'left'
    new_boundary = 'nsr'
  []
  [top]
    type = CartesianMeshGenerator
    dim = 2
    dx = '6 1 3'
    dy = '5 1 3'
    ix = '12 2 6'
    iy = '10 2 6'
    subdomain_id = '0 0 0
                    0 1 0
                    0 0 0'
  []
  [top_trans]
    type = TransformGenerator
    input = top
    transform = TRANSLATE
    vector_value = '0.0 1.0 0.0'
  []
  [top_rename]
    type = RenameBoundaryGenerator
    input = top_trans
    old_boundary = 'bottom'
    new_boundary = 'nsr'
  []
  [all]
    type = StitchedMeshGenerator
    inputs = 'top_rename bottom_rename'
    clear_stitched_boundary_ids = true
    stitch_boundaries_pairs = 'nsr top'
    merge_boundaries_with_same_name = false
  []
  [merge_1]
    type = RenameBoundaryGenerator
    input = all
    old_boundary = '1 5 0 7'
    new_boundary = 'right right nsr nsr'
  []
  [merge_2]
    type = RenameBoundaryGenerator
    input = merge_1
    old_boundary = 'right top left bottom'
    new_boundary = 'vacuum vacuum vacuum vacuum'
  []
  uniform_refine = 0
[]
