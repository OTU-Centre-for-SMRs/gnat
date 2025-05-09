# Simple 2D Fixed-Source Problems

In this tutorial, you will learn:

1. How to model 2D shielding problems in Gnat;
2. How to provide multi-group material properties;
2. How to specify point and surface sources.

To access this tutorial,

```bash
cd gnat/tutorials/2D_fixed_source
```

## Geometry and Material Properties

This problem is one where an isotropic neutron point source is placed in the
middle of a $10\text{ cm }\times 10\text{ cm }$ 2D domain that both absorbs and scatters neutrons isotropically.
The point source emits at a constant rate of $1\times 10^{3}\text{ s}^{-1}$ into a single energy group. Vacuum
boundaries surround the domain on all four sides and the problem is at
steady-state. The material has $\Sigma_{t} = 2.0$ cm$^{-1}$ and $\Sigma_{s} = 1.0$ cm$^{-1}$.

## The Neutronics Input File

Like any problem which is solved with the finite element method, the first step
in using Gnat for neutron transport is to declare a mesh. For a simple case
like this one the built-in mesh generators in Moose can be used:

!listing /tutorials/2D_fixed_source/2D_point_source.i
  block=Mesh

The [CartesianMeshGenerator](https://mooseframework.inl.gov/source/meshgenerators/CartesianMeshGenerator.html)
object generates simple cartesian meshes in 1D, 2D or 3D. `dim = 2` is set to let
the mesh generator know that a 2D cartesian mesh is required. `dx = 10` and
`dy = 10` set the size of the mesh to 10 cm in both x and y. Finally, `ix = 100`
and `iy = 100` indicate that the mesh should be subdivided into 100 elements
along both the x and y axis. The result is shown in [!ref](2D_problem_mesh):

!media media/2D_fixed_source/mesh.png id=2D_problem_mesh caption=2D cartesian mesh for the 2D point source problem.
  style=width:40%;margin-left:auto;margin-right:auto;halign:center

The next step in solving radiation transport problems is to declare a `TransportSystem` block. This is short-cut
syntax that sets up all of the kernels and variables necessary for transport calculations along many directions
for many energy groups, which minimizes the length of your Gnat input file.

!listing /tutorials/2D_fixed_source/2D_point_source.i
  block=TransportSystems

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
1000 n/s emitting into the single energy group.

Next, the material properties of the domain must be specified. Gnat does this
through a system implemented for transport problems known as 'TransportMaterials':

!listing /tutorials/2D_fixed_source/2D_point_source.i
  block=TransportMaterials

This declares a constant material which provides enough parameters for a single
neutron group (`num_groups = 1`). The material provides an isotropic scattering
cross section (`anisotropy = 0`). The material has an total cross section
of 2 (`group_total = 2.0`) and a scattering cross section of 1
(`group_scattering = 1.0`). `transport_system = Neutron` links this material to
the `Neutron` transport system previously declared.

The remainder of the input file describes how the problem should be solved:

!listing /tutorials/2D_fixed_source/2D_point_source.i
  block=Executioner

The `Executioner` block declares the solver parameters and solve type. In this case,
`type = Steady` and the Preconditioned Jacobian-Free Newton Krylov
(`solve_type = PJFNK`) is used. It is recommended that the user employ the PJFNK
solver for continuous finite element methods when solving radiation transport
problems with the SAAF-CFEM scheme.

The simulation can be executed with the following shell command:

```language=bash
mpiexec -np 2 gnat-opt -i ./2D_point_source.i --n-threads=2
```

The results of this simple case can be seen below in [!ref](point_output):

!media media/2D_fixed_source/1_group_point_source.png id=point_output caption=Scalar flux from the 2D point source problem.
  style=width:40%;margin-left:auto;margin-right:auto;halign:center

## Adding a Second Energy Group

## Swapping the Point Source for a Surface Source
