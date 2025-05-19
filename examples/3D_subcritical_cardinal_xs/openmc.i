[Mesh]
  [file]
    type = FileMeshGenerator
    file = mesh_in.e
  []
[]

[AuxVariables]
  [Scatter_Ratio_G1]
    type = MooseVariable
    family = MONOMIAL
    order = CONSTANT

    block = 'fuel cladding graphite_core graphite_reflector stainless_steel
             fuel_trimmer_tri graphite_core_trimmer_tri graphite_reflector_trimmer_tri
             stainless_steel_trimmer_tri'
  []
  [Scatter_Ratio_G2]
    type = MooseVariable
    family = MONOMIAL
    order = CONSTANT

    block = 'fuel cladding graphite_core graphite_reflector stainless_steel
             fuel_trimmer_tri graphite_core_trimmer_tri graphite_reflector_trimmer_tri
             stainless_steel_trimmer_tri'
  []
  [Scatter_Ratio_G3]
    type = MooseVariable
    family = MONOMIAL
    order = CONSTANT

    block = 'fuel cladding graphite_core graphite_reflector stainless_steel
             fuel_trimmer_tri graphite_core_trimmer_tri graphite_reflector_trimmer_tri
             stainless_steel_trimmer_tri'
  []
[]

[AuxKernels]
  [Comp_Scatter_Ratio_G1]
    type = ParsedAux
    variable = Scatter_Ratio_G1
    coupled_variables = 'total_xs_g1 scatter_xs_g1_gp1_l0 scatter_xs_g1_gp2_l0 scatter_xs_g1_gp3_l0'
    expression = '(scatter_xs_g1_gp1_l0 + scatter_xs_g1_gp2_l0 + scatter_xs_g1_gp3_l0) / total_xs_g1'

    block = 'fuel cladding graphite_core graphite_reflector stainless_steel
             fuel_trimmer_tri graphite_core_trimmer_tri graphite_reflector_trimmer_tri
             stainless_steel_trimmer_tri'
  []
  [Comp_Scatter_Ratio_G2]
    type = ParsedAux
    variable = Scatter_Ratio_G2
    coupled_variables = 'total_xs_g2 scatter_xs_g2_gp1_l0 scatter_xs_g2_gp2_l0 scatter_xs_g2_gp3_l0'
    expression = '(scatter_xs_g2_gp1_l0 + scatter_xs_g2_gp2_l0 + scatter_xs_g2_gp3_l0) / total_xs_g2'

    block = 'fuel cladding graphite_core graphite_reflector stainless_steel
             fuel_trimmer_tri graphite_core_trimmer_tri graphite_reflector_trimmer_tri
             stainless_steel_trimmer_tri'
  []
  [Comp_Scatter_Ratio_G3]
    type = ParsedAux
    variable = Scatter_Ratio_G3
    coupled_variables = 'total_xs_g3 scatter_xs_g3_gp1_l0 scatter_xs_g3_gp2_l0 scatter_xs_g3_gp3_l0'
    expression = '(scatter_xs_g3_gp1_l0 + scatter_xs_g3_gp2_l0 + scatter_xs_g3_gp3_l0) / total_xs_g3'

    block = 'fuel cladding graphite_core graphite_reflector stainless_steel
             fuel_trimmer_tri graphite_core_trimmer_tri graphite_reflector_trimmer_tri
             stainless_steel_trimmer_tri'
  []
[]

[Problem]
  type = OpenMCCellAverageProblem

  particles = 1000
  batches = 600
  inactive_batches = 50

  verbose = false
  power = 1.0
  source_rate_normalization = 'kappa_fission'
  cell_level = 1

  # We aren't tallying on the regions that contain air, which results in a few missed hits.
  check_tally_sum = false

  [MGXS]
    tally_type = cell
    particle = neutron
    energy_boundaries = '14900000.0 111000.0 2.38 0.0'
    estimator = 'analog'
    block = 'fuel cladding graphite_core graphite_reflector stainless_steel
             fuel_trimmer_tri graphite_core_trimmer_tri graphite_reflector_trimmer_tri
             stainless_steel_trimmer_tri'

    add_scattering = true
    transport_correction = false

    add_fission = true
    add_fission_heating = true

    hide_tally_vars = false
  []
[]

[Executioner]
  type = Steady
[]

[Outputs]
  exodus = true
  csv = true
  execute_on = TIMESTEP_END
[]
