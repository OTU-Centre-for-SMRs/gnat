# Mesh is in cm, need to convert all units to cm.
[Mesh]
  [Domain]
    type = FileMeshGenerator
    file = Radiation_Mesh.e
  []
[]

[TransportSystems]
  [Photon]
    scheme = saaf_cfem
    particle_type = photon

    flux_moment_names = 'Flux_Moment'
    family = LAGRANGE
    order = FIRST

    num_groups = 1
    max_anisotropy = 0
    n_azimuthal = 3
    n_polar = 3

    vacuum_boundaries = vacuum

    field_source_blocks = air
    field_source_moments = 'photon_source_g_0'
    field_source_anisotropies = '0'
    field_source_scaling = '1e5'
    scale_sources = true
  []
[]

[TransportMaterials]
  [Air]
    type = ConstantTransportMaterial
    transport_system = Photon
    group_total = 3.29e-5
    group_scattering = '4.6e-6'
    anisotropy = 0
  []
[]

[AuxVariables]
  [photon_source_g_0]
    type = MooseVariable
    family = MONOMIAL
    order = CONSTANT
  []
[]

[Transfers]
  [G0]
    type = MultiAppGeneralFieldShapeEvaluationTransfer
    from_multi_app = Tracer
    variable = photon_source_g_0
    source_variable = photon_source_g_0
  []
[]

[MultiApps]
  [Tracer]
    type = TransientMultiApp
    app_type = GnatApp
    execute_on = TIMESTEP_BEGIN
    input_files = tracer_transient.i
  []
[]

[Executioner]
  type = Transient
  solve_type = PJFNK

  petsc_options_iname = '-ksp_gmres_restart'
  petsc_options_value = ' 600'

  [TimeStepper]
    type = ConstantDT
    dt = 5.0
  []
  end_time = 60

  l_max_its = 50
  nl_abs_tol = 1e-12

  automatic_scaling = true
  off_diagonals_in_auto_scaling = true
  compute_scaling_once = false
[]

[Outputs]
  exodus = true
[]
