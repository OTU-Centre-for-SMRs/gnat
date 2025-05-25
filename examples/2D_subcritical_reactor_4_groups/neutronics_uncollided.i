# An input deck which compute the uncollided flux in the assembly with the SASF approach.
# The SASF method in 2D assumes a plane-wave point source which is not physically realistic for most
# particle sources.

[Mesh]
  [NeutronicsDomain]
    type = FileMeshGenerator
    file = 2D_quarter_core_nsr.e
  []
  uniform_refine = 0
[]

[UncollidedFlux]
  [Neutron]
    uncollided_flux_treatment = sasf
    num_groups = 4
    max_anisotropy = 0

    point_source_locations = '0.0 0.0 0.0'
    point_source_moments = '1.5'
    point_source_anisotropies = '0'

    sasf_near_source_boundary = 'near_source'
    sasf_vacuum_boundaries = 'vacuum'
    sasf_near_source_cross_sections = '9.683916522303286e-05 0.00023453400588691886 0.00043185580016127494 0.0005108576781703109'

    is_conservative_transfer_src = false
  []
[]

[TransportMaterials]
  [Air]
    type = FileTransportMaterial
    transport_system = Neutron
    file_name = '../../data/mgxs/sc_4g_xs_macro.xml'
    source_material_id = '5' #
    block = 'air_gap'
  []

  [Cladding]
    type = FileTransportMaterial
    transport_system = Neutron
    file_name = '../../data/mgxs/sc_4g_xs_macro.xml'
    source_material_id = '7' #
    block = 'cladding'
  []

  [Fuel]
    type = FileTransportMaterial
    transport_system = Neutron
    file_name = '../../data/mgxs/sc_4g_xs_macro.xml'
    source_material_id = '6' #
    block = 'fuel'
  []

  [Graphite]
    type = FileTransportMaterial
    transport_system = Neutron
    file_name = '../../data/mgxs/sc_4g_xs_macro.xml'
    source_material_id = '8' #
    block = 'core_graphite reflector_graphite'
  []

  [Control]
    type = FileTransportMaterial
    transport_system = Neutron
    file_name = '../../data/mgxs/sc_4g_xs_macro.xml'
    source_material_id = '8' #
    block = 'control_or_graphite'
  []
[]

[Executioner]
  type = Steady
  solve_type = PJFNK
  petsc_options_iname = '-ksp_gmres_restart'
  petsc_options_value = ' 600'

  l_max_its = 50
  nl_abs_tol = 1e-12

  automatic_scaling = true
  off_diagonals_in_auto_scaling = true
  compute_scaling_once = false
[]
