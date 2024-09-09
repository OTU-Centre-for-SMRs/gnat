# A Simple Neutron Transport Problem

This problem is one where an isotropic neutron point source is placed in the
middle of a 2D domain that both absorbs and scatters neutrons isotropically.
The point source emits at a constant rate into a single energy group. Vacuum
boundaries surround the domain on all four sides and the problem is at
steady-state.

Like any problem which is solved with the finite element method, the first step
in using Gnat for neutron transport is to declare a mesh. For a simple case
like this one the built-in mesh generators in Moose can be used:

```language=moose
[Mesh]
  [domain]
    type = CartesianMeshGenerator
    dim = 2
    dx = 10
    dy = 10
    ix = 100
    iy = 100
  []
[]
```

[/examples/2D_scattering_cgfem/test.i]

The `CartesianMeshGenerator` object generates simple cartesian meshes in 1D, 2D
or 3D. `dim = 2` is set to let the mesh generator know that a 2D cartesian mesh
is required. `dx = 10` and `dy = 10` set the size of the mesh to 10 cm in both x
and y. Finally, `ix = 100` and `iy = 100` indicate that the mesh should be
subdivided into 100 elements along both the x and y axis. The result is shown
in [!ref](example_1_mesh):

!media media/example_1/example_1_mesh.png id=example_1_mesh caption=2D cartesian mesh for Example 1.
  style=width:40%;margin-left:auto;margin-right:auto;halign:center

The next step in solving radiation transport problems is to declare a `TransportSystem` block. These are custom Moose
objects added by Gnat to make setting up radiation transport simulations as
pain-free as possible.

```language=moose
[TransportSystems]
  [Neutron]
    num_groups = 1
    max_anisotropy = 0
    scheme = saaf_cfem
    particle_type = neutron

    order = FIRST
    family = LAGRANGE

    n_azimuthal = 2
    n_polar = 2

    vacuum_boundaries = 'left right top bottom'

    point_source_locations = '5.0 5.0 0.0'
    point_source_moments = '1e3'
    point_source_anisotropies = '0'
    scale_sources = true
  []
[]
```

[/examples/2D_scattering_cgfem/test.i]

To start, `num_groups = 1` specifies that the transport simulation should
only consider a single neutron group although more can be specified at the cost
of increased computational complexity. `max_anisotropy = 0` tells the transport
simulation to only evaluate the 0th-order moment of the angular flux, which is
more than enough for this case as scattering is isotropic. `scheme = saaf_cfem`
states that the self-adjoint angular flux stabilization scheme is requested.
In accordance with that choice, continuous basis functions in the form of
first-order (`order = FIRST`) Lagrange (`family = LAGRANGE`) functions are used.
The number of angular directions are specified with `n_azimuthal = 2` and
`n_polar = 2`, and vacuum boundaries are specified for all sides with
`vacuum_boundaries = 'left right top bottom'`. Finally, the point source is
specified with a location at (5.0, 5.0, 5.0) with an emission intensity of
1000 n/s emitting into group 1.

Next, the material properties of the domain must be specified. Gnat does this
through a system implemented for transport problems known as 'TransportMaterials':

```language=moose
[TransportMaterials]
  [Domain]
    type = ConstantTransportMaterial
    transport_system = Neutron
    anisotropy = 0
    group_total = 2.0
    group_scattering = 1.0
  []
[]
```

[/examples/2D_scattering_cgfem/test.i]

This declares a constant material which provides enough parameters for a single
neutron group (`num_groups = 1`). The material provides an isotropic scattering
cross section (`anisotropy = 0`). The material has an total cross section
of 2 (`group_total = 2.0`) and a scattering cross section of 1
(`group_scattering = 1.0`). `transport_system = Neutron` links this material to
the `Neutron` transport system previously declared.

The remainder of the input file describes how the problem should be solved:

```language=moose
[Executioner]
  type = Steady
  solve_type = PJFNK
  petsc_options_iname = '-pc_type -pc_hypre_type -ksp_gmres_restart'
  petsc_options_value = 'hypre boomeramg 10'
  l_max_its = 50
  nl_rel_tol = 1e-12
[]
```

[/examples/2D_scattering_cgfem/test.i]

The `Executioner` block declares the solver parameters and solve type. In this case,
`type = Steady` and the Preconditioned Jacobian-Free Newton Krylov
(`solve_type = PJFNK`) is used. It is recommended that the user employ the PJFNK
solver for continuous finite element methods when solving radiation transport
problems dominated by large scattering ratios.

The simulation can be executed with the following shell command:

```language=bash
./gnat-opt -i ./examples/2D_scattering_cgfem/test.i
```

The results of this simple case can be seen below in [!ref](example_1_output):

!media media/example_1/example_1_output.png id=example_1_output caption=Example 1 Results.
  style=width:40%;margin-left:auto;margin-right:auto;halign:center
