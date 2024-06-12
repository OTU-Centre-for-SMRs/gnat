[Mesh]
  type = GeneratedMesh
  dim = 2
  nx = 8
  ny = 8
[]

[Variables]
  [u]
  []
[]

[Kernels]
  [Diffusion]
    type = ADMassFractionNuclideDiffusion
    variable = u
    u = 1.0
    v = 1.0
    density = 1.0
  []

  [Advection]
    type = ADMassFractionNuclideAdvection
    variable = u
    u = 1.0
    v = 1.0
    density = 1.0
  []

  [Decay]
    type = ADMassFractionNuclideDecaySink
    variable = u
    u = 1.0
    v = 1.0
    decay_const = 1.0
    density = 1.0
  []

  [Force]
    type = ADIsotopeForcing
    variable = u
    u = 1.0
    v = 1.0
    forcing = force
    density = 1.0
  []
[]

[BCs]
  [all]
    type = FunctionDirichletBC
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
  type = Steady
  solve_type = 'PJFNK'
  petsc_options_iname = '-pc_type -pc_hypre_type'
  petsc_options_value = 'hypre boomeramg'
[]

[Outputs]
  exodus = true
  csv = true
[]
