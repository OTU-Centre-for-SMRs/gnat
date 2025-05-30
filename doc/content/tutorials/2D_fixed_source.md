# Simple Fixed-Source Problems

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
middle of a $10\text{ cm }\times 10\text{ cm }$ 2D domain that scatters neutrons isotropically and absorbs neutrons.
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
1000 n/s emitting into the single energy group. In addition to specifying the
point source, we also set `scale_sources = true` which divides all source moments
by the largest moment value. This rescales the fixed source problem such that
the largest source moment is unity, improving the stability of the finite element
solve.

Next, the material properties of the domain must be specified. Gnat does this
through a system implemented for transport problems known as `TransportMaterials`:

!listing /tutorials/2D_fixed_source/2D_point_source.i
  block=TransportMaterials

This declares a constant material which provides enough parameters for a single
neutron group. The material provides an isotropic scattering
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

We can see the signature disadvantage of the discrete ordiantes method: ray effects. These
occur because of the limited angular resolution (four directions per quadrant of the 2D unit sphere)
of the chosen discretization compared to the small size of the fixed source. We can attempt to mitigate
these ray effects by adding extra directions by specifying `n_azimuthal = 10` and `n_polar = 10`, which
results in [point_output_ref].

!media media/2D_fixed_source/1_group_point_source_ref.png id=point_output_ref caption=Refined scalar flux from the 2D point source problem.
  style=width:40%;margin-left:auto;margin-right:auto;halign:center

However, the increased angular fidelity results in a large increase in runtime with diminishing returns.
Instead of increasing the number of directions one could instead use an uncollided flux technique,
which is discussed in the [uncollided flux tutorial](uncollided_flux.md).

## Adding a Second Energy Group

In this section, we modify the problem to add a second energy group. The modifications
made to the transport system can be found below.

!listing /tutorials/2D_fixed_source/2D_point_source_2_grp.i
  block=TransportSystems

We specify that two energy groups should be added by setting `num_groups = 2`. We then add
a second source moment for the point source with `point_source_moments = '1e3 1.0'`. The
`TransportMaterials` block then needs to be modified with a new set of group properties.

!listing /tutorials/2D_fixed_source/2D_point_source_2_grp.i
  block=TransportMaterials

Here, we modify the total cross section such that it's representative of two-group transport in
a highly scattering material (such as water). The scattering cross section becomes a scattering matrix
with zero thermal upscattering. This modified input deck can be run with:

```language=bash
mpiexec -np 2 gnat-opt -i ./2D_point_source_2_grp.i --n-threads=2
```

The results for group 1 (fast) can be found in [2g_point_g1], while the results for group 2 (thermal) can be
found in [2g_point_g2].

!row! style=display:inline-flex;
!col!

!media media/2D_fixed_source/2_group_point_source_g1.png id=2g_point_g1 caption=Group 1 scalar flux from the 2D point source problem.
  style=width:100%;margin-left:auto;margin-right:auto;halign:center

!col-end!

!col!

!media media/2D_fixed_source/2_group_point_source_g2.png id=2g_point_g2 caption=Group 2 scalar flux from the 2D point source problem.
  style=width:100%;margin-left:auto;margin-right:auto;halign:center

!col-end!

!row-end!

Similar to [point_output], we can see ray effects in the fast energy group due to the low within-group scattering
cross section relative to the total cross section. The effect of fast-to-thermal scattering and the large
within-group thermal scattering cross section can be seen in [2g_point_g2], where there are no ray effect
even though the point source emits into group two.

## Swapping the Point Source for a Surface Source

In this section, we remove the point source and add an isotropic surface source on the top boundary of the domain.
This requires a few minor modifications to the `TransportSystem`.

!listing /tutorials/2D_fixed_source/2D_surface_source_2_grp.i
  block=TransportSystems

`point_source_locations` is replaced with `source_boundaries`, which contains the side sets on which we wish to apply the
surface source. `point_source_moments` and `point_source_anisotropies` are then replaced with their surface source equivalents:
`boundary_source_moments` and `boundary_source_anisotropy`. This modified input deck can be run with:

```language=bash
mpiexec -np 2 gnat-opt -i ./2D_surface_source_2_grp.i --n-threads=2
```

The results for group 1 (fast) can be found in [2g_surf_g1], while the results for group 2 (thermal) can be
found in [2g_surf_g2].

!row! style=display:inline-flex;
!col!

!media media/2D_fixed_source/2_group_surface_source_g1.png id=2g_surf_g1 caption=Group 1 scalar flux from the 2D surface source problem.
  style=width:100%;margin-left:auto;margin-right:auto;halign:center

!col-end!

!col!

!media media/2D_fixed_source/2_group_surface_source_g2.png id=2g_surf_g2 caption=Group 2 scalar flux from the 2D surface source problem.
  style=width:100%;margin-left:auto;margin-right:auto;halign:center

!col-end!

!row-end!
