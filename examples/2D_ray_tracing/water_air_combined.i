# A simple test case with a purely absorbing medium and a point source in the
# middle of the domain.

[Mesh]
  [domain]
    type = CartesianMeshGenerator
    dim = 2
    dx = '5.0 10.0'
    dy = '10.0 5.0'
    ix = '5 10'
    iy = '10 5'
    subdomain_id = '
      1 1
      2 1'
  []
  uniform_refine = 1
[]

[TransportSystems]
  [Neutron]
    num_groups = 2
    max_anisotropy = 1
    scheme = saaf_cfem
    particle_type = neutron

    order = FIRST
    family = LAGRANGE

    n_azimuthal = 5
    n_polar = 5

    vacuum_boundaries = 'left right top bottom'

    uncollided_from_multi_app = 'Uncollided_Neutronics'
  []
[]

[TransportMaterials]
  [Water]
    type = FileNeutronicsMaterial
    transport_system = Neutron
    file_name = 'macro_xs_water.xml'
    source_material_id = '1'
    block = '1'
  []
  [Air]
    type = FileNeutronicsMaterial
    transport_system = Neutron
    file_name = 'macro_xs_air.xml'
    source_material_id = '5'
    block = '2'
  []
[]

[Executioner]
  type = Steady
  solve_type = PJFNK
  #petsc_options_iname = '-pc_type -pc_hypre_type -ksp_gmres_restart'
  #petsc_options_value = ' hypre    boomeramg      100'
  l_max_its = 50
  nl_rel_tol = 1e-12

  automatic_scaling = true
  off_diagonals_in_auto_scaling = true
  compute_scaling_once = false
[]

[Outputs]
  exodus = true
[]

[MultiApps]
  [Uncollided_Neutronics]
    type = FullSolveMultiApp
    app_type = GnatApp
    input_files = 'water_air_unc.i'
    execute_on = INITIAL
  []
[]
