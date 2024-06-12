[Mesh]
  type = GeneratedMesh
  dim = 2
  nx = 64
  ny = 64
[]

[Variables]
  [u]
    type = MooseVariableFVReal
    fv = true
  []
[]

[FVKernels]
  [Time]
    type = ADFVMassFractionTimeDerivative
    variable = u
    density = 1.0
  []

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
    advected_interp_method = vanLeer
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
    value = 'x*y*alpha*t^3 + 3*x*y*t^2 + x*t^3*vec_y + y*t^3*vec_x'
    vars = 'vec_x alpha vec_y'
    vals = '1.0 1.0 1.0'
  []
  [exact]
    type = ParsedFunction
    value = 'x*y*t^3'
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
[]

[Executioner]
  type = Transient
  solve_type = NEWTON
  petsc_options_iname = '-pc_type -pc_factor_shift_type'
  petsc_options_value = ' lu       NONZERO'

  dt = 1
  end_time = 3

  line_search = 'none'

  automatic_scaling = true
  off_diagonals_in_auto_scaling = true
  compute_scaling_once = false
[]

[Outputs]
  exodus = true
  csv = true
[]
