# A simple test case with a purely absorbing medium and a point source in the
# middle of the domain.

[Mesh]
  [gen]
    type = GeneratedMeshGenerator
    dim = 2
    xmin = 0
    xmax = 10
    ymin = 0
    ymax = 10
    nx = 101
    ny = 101
  []
[]

[TransportSystems]
  [Neutron]
    num_groups = 2
    max_anisotropy = 0
    scheme = saaf_cfem
    particle_type = neutron

    order = FIRST
    family = LAGRANGE

    n_azimuthal = 10
    n_polar = 10

    vacuum_boundaries = 'left right top bottom'

    point_source_locations = '5.0 5.0 0.0'
    point_source_moments = '1.0 0.0'
    point_source_anisotropies = '0'
    scale_sources = true
  []
[]

[TransportMaterials]
  [Domain]
    type = FileNeutronicsMaterial
    transport_system = Neutron
    file_name = 'macro_xs.xml'
    source_material_id = '5'
  []
[]

[Executioner]
  type = Steady
  solve_type = PJFNK

  nl_abs_tol = 1e-12
  l_max_its = 50

  automatic_scaling = true
  off_diagonals_in_auto_scaling = true
  compute_scaling_once = false
[]

[Outputs]
  exodus = true
[]
