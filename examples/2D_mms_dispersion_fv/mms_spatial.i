[Mesh]
  type = GeneratedMesh
  dim = 2
  nx = 8
  ny = 8
[]

[Variables]
  [u]
    type = MooseVariableFVReal
    fv = true
  []
[]

[FVKernels]
  [Diffusion]
    type = ADFVMassFractionNuclideDiffusion
    variable = u
    density = 1.0
  []

  [Advection]
    type = ADFVMassFractionFunctorAdvection
    variable = u
    density = 1.0
    u = 1.0
    v = 1.0
    advected_interp_method = upwind
    velocity_interp_method = average
  []

  [Decay]
    type = ADFVMassFractionNuclideDecaySink
    variable = u
    density = 1.0
    decay_const = 1.0
  []

  [Force]
    type = FVBodyForce
    variable = u
    function = force
  []
[]

[FVBCs]
  [All]
    type = FVFunctionDirichletBC
    variable = u
    function = exact
    boundary = 'left right top bottom'
  []
[]

[Materials]
  [Mat1]
    type = SUPGFunctorTestMaterial
    var_name = u
    diff = 1.0
    vel = '1.0 1.0 0.0'
  []
[]

[Functions]
  [force]
    type = ParsedFunction
    value = '-x^2*D*exp(x*y) + x*vec_y*exp(x*y) - y^2*D*exp(x*y) + y*vec_x*exp(x*y) + alpha*exp(x*y)'
    vars = 'alpha vec_y D vec_x'
    vals = '1.0 1.0 1.0 1.0'
  []
  [exact]
    type = ParsedFunction
    value = 'exp(x*y)'
  []
  [peclet]
    type = ParsedFunction
    value = 'h * sqrt(2.0) / 1.0'
    symbol_names = 'h'
    symbol_values = 'h'
  []
[]

[Postprocessors]
  [error]
    type = ElementL2Error
    function = exact
    variable = u
  []
  [h]
    type = AverageElementSize
  []
  [PecletNumber]
    type = FunctionValuePostprocessor
    function = peclet
    indirect_dependencies = 'h'
  []
[]

[Executioner]
  type = Steady
  solve_type = NEWTON
  petsc_options_iname = '-pc_type -pc_factor_shift_type'
  petsc_options_value = ' lu       NONZERO'

  line_search = 'none'

  automatic_scaling = true
  off_diagonals_in_auto_scaling = true
  compute_scaling_once = false
[]

[Outputs]
  exodus = true
  csv = true
[]
