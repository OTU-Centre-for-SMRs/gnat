[Mesh]
  type = GeneratedMesh
  dim = 2
  nx = 64
  ny = 64
[]

[Variables]
  [u]
  []
[]

[Kernels]
  [Time]
    type = ADMassFractionNuclideTimeDerivative
    variable = u
    density = 1.0
    u = 1.0
    v = 1.0
  []

  [Diffusion]
    type = ADMassFractionNuclideDiffusion
    variable = u
    density = 1.0
    u = 1.0
    v = 1.0
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
    decay_const = 1.0
    density = 1.0
    u = 1.0
    v = 1.0
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
    value = 'x*y*alpha*t^3 + 3*x*y*t^2 + x*t^3*vec_y + y*t^3*vec_x'
    vars = 'alpha vec_y vec_x'
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
  dt = 1
  end_time = 3
  solve_type = 'PJFNK'
  petsc_options_iname = '-pc_type -pc_hypre_type'
  petsc_options_value = 'hypre boomeramg'
[]

[Outputs]
  exodus = true
  csv = true
[]
