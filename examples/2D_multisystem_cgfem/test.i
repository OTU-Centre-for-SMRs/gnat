# A simple test case with a purely absorbing medium and a point source in the
# middle of the domain.

[Mesh]
  [domain]
    type = CartesianMeshGenerator
    dim = 2
    dx = '10.0'
    dy = '10.0'
    ix = '101'
    iy = '101'
    subdomain_id = '
      1'
  []
[]

[TransportSystems]
  [NeutronSN]
    scheme = saaf_cfem
    particle_type = neutron
    num_groups = 2
    max_anisotropy = 0
    angular_flux_names = neutron_angular_flux
    flux_moment_names = neutron_flux_moment

    order = FIRST
    family = LAGRANGE

    n_azimuthal = 3
    n_polar = 3

    vacuum_boundaries = 'left right top bottom'

    point_source_locations = '5.0 5.0 0.0'
    point_source_moments = '1.0 0.0'
    point_source_anisotropies = '0'
    scale_sources = true
  []

  [PhotonSN]
    scheme = saaf_cfem
    particle_type = Photon
    num_groups = 2
    max_anisotropy = 0
    angular_flux_names = photon_angular_flux
    flux_moment_names = photon_flux_moment

    order = FIRST
    family = LAGRANGE

    n_azimuthal = 3
    n_polar = 3

    vacuum_boundaries = 'left right top bottom'

    point_source_locations = '5.0 5.0 0.0'
    point_source_moments = '0.0 1.0'
    point_source_anisotropies = '0'
    scale_sources = true
  []
[]

[TransportMaterials]
  [Water1]
    type = FileTransportMaterial
    transport_system = NeutronSN
    file_name = '../../data/mgxs/water_2g_xs_macro.xml'
    source_material_id = '1'
    block = '1'
  []

  [Water2]
    type = FileTransportMaterial
    transport_system = PhotonSN
    file_name = '../../data/mgxs/water_2g_xs_macro.xml'
    source_material_id = '1'
    block = '1'
  []
[]

[Outputs]
  exodus = true
[]

[Executioner]
  type = Steady
  solve_type = PJFNK
  petsc_options_iname = '-pc_type -pc_hypre_type -ksp_gmres_restart'
  petsc_options_value = ' hypre    boomeramg      600'

  l_max_its = 50
  nl_rel_tol = 1e-12
[]
