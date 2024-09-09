[Mesh]
  [bottom_1]
    type = CartesianMeshGenerator
    dim = 3
    dx = '50'
    dy = '10'
    dz = '10'
    ix = '5'
    iy = '1'
    iz = '1'
    subdomain_id = '0'
  []
  [bottom_1_trans]
    type = TransformGenerator
    input = bottom_1
    transform = TRANSLATE
    vector_value = '10.0 0.0 0.0'
  []
  [bottom_1_rename]
    type = RenameBoundaryGenerator
    input = bottom_1_trans
    old_boundary = 'left'
    new_boundary = 'nsr'
  []
  [bottom_2]
    type = CartesianMeshGenerator
    dim = 3
    dx = '10 20 10 20'
    dy = '40 10 40'
    dz = '10'
    ix = '1 2 1 2'
    iy = '4 1 4'
    iz = '1'
    subdomain_id = '1 0 0 0
                    1 1 1 0
                    0 0 0 0'
  []
  [bottom_2_trans]
    type = TransformGenerator
    input = bottom_2
    transform = TRANSLATE
    vector_value = '0.0 10.0 0.0'
  []
  [bottom_2_rename]
    type = RenameBoundaryGenerator
    input = bottom_2_trans
    old_boundary = 'bottom'
    new_boundary = 'nsr'
  []
  [stitch_bottom]
    type = StitchedMeshGenerator
    inputs = 'bottom_2_rename bottom_1_rename'
    clear_stitched_boundary_ids = true
    stitch_boundaries_pairs = 'nsr top'
    merge_boundaries_with_same_name = false
  []
  [bottom_rename_all]
    type = RenameBoundaryGenerator
    input = stitch_bottom
    old_boundary = '1 10 0 6 2 8 3 4 7 5 11'
    new_boundary = 'nsr nsr vacuum vacuum vacuum vacuum vacuum vacuum vacuum front front'
  []
  [top]
    type = CartesianMeshGenerator
    dim = 3
    dx = '10 20 10 20'
    dy = '10 40 10 40'
    dz = '20 10 20'
    ix = '1 2 1 2'
    iy = '1 4 1 4'
    iz = '2 1 2'
    subdomain_id = '0 0 0 0
                    0 0 0 0
                    0 0 1 0
                    0 0 0 0

                    0 0 0 0
                    0 0 0 0
                    0 0 1 0
                    0 0 1 0

                    0 0 0 0
                    0 0 0 0
                    0 0 0 0
                    0 0 0 0'
  []
  [top_trans]
    type = TransformGenerator
    input = top
    transform = TRANSLATE
    vector_value = '0.0 0.0 10.0'
  []
  [top_rename_all]
    type = RenameBoundaryGenerator
    input = top_trans
    old_boundary = 'left right front back top bottom'
    new_boundary = 'vacuum vacuum vacuum nsr vacuum vacuum'
  []
  [stitch_all]
    type = StitchedMeshGenerator
    inputs = 'top_rename_all bottom_rename_all'
    clear_stitched_boundary_ids = true
    stitch_boundaries_pairs = 'nsr front'
    merge_boundaries_with_same_name = false
  []
  [rename_all]
    type = RenameBoundaryGenerator
    input = stitch_all
    old_boundary = '0 6 4 7'
    new_boundary = 'nsr nsr vacuum vacuum'
  []
[]
