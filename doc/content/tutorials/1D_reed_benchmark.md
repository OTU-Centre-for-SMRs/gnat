# The Reed Benchmark Problem

In this tutorial, you will learn:

1. The basics of setting up a neutron transport problems in Gnat;
2. How to model the 1D Reed benchmark problem.

To access this tutorial,

```bash
cd gnat/tutorials/1D_reed_sn
```

## Geometry and Material Properties

The Reed problem is a classic benchmark proposed in [!cite](reed) which is used to test the accuracy of space-angle
discretization schemes, and serves as the "hello world" in the radiation transport community. It consists of a
1D slab with five distinct material regions, which can be found below in [table1].

!table id=table1 caption=Material properties for the 1D Reed problem.
| ID | Region | $\Sigma_{t}$ ($\text{cm}^{-1}$) | $\Sigma_{s}$ ($\text{cm}^{-1}$) | $q$ ($\text{s}^{-1}$) |
| :- | - | - | - | - |
| 1 | $0\leq x < 2.0$ | 50.0 | 0.0 | 50.0 |
| 2 | $2\leq x < 3.0$ | 5.0 | 0.0 | 50.0 |
| 3 | $3\leq x < 5.0$ | 0.0 | 0.0 | 50.0 |
| 4 | $5\leq x < 6.0$ | 1.0 | 0.9 | 1.0 |
| 5 | $6\leq x < 8.0$ | 1.0 | 0.9 | 0.0 |

The slab is arranged such that the left boundary is a reflective (symmetry) boundary,
and the right boundary is a vacuum boundary condition.

## The Neutronics Input File

Like any problem which is solved with the finite element method, the first step
when using Gnat for neutron transport is to declare a mesh. The Reed benchmark is
incredibly simple, so we can use MOOSE's built-in [CartesianMeshGenerator](https://mooseframework.inl.gov/source/meshgenerators/CartesianMeshGenerator.html):

!listing /tutorials/1D_reed_sn/1D_reed.i
  block=Mesh

Here, we specify the mesh should be 1D with `dim = 1`. This is
followed by the unique spatial regions in `dx`. Afterwards, the number of subdivisions per
region are specified in `ix` - initially set such that that each element is $0.5$ cm long.
Finally, we set the unique material regions in `blocks`. This mesh generator automatically
generates side sets called `left` and `right`, which we can use to apply boundary conditions.

The next step is to set up the radiation transport problem. If we were to do this with
MOOSE's default input syntax, we would need to add variables and kernels for each direction and
each energy group. This would be incredibly tedious, so Gnat provides input syntax to automate
this process with the [TransportSystems block](TransportAction.md):

!listing /tutorials/1D_reed_sn/1D_reed.i
  block=TransportSystems

Here, we create a "transport system" named `Neutron`. We specify that it should be mono-energetic
with `num_groups = 1`, and that it should use the SAAF-CFEM-SN scheme for discretizing the transport
equation (`scheme = saaf_cfem`). We then specify that `particle_type = neutron`, indicating that we
want to add kernels relevant to neutrons. For the purposes of the Reed benchmark problem, one could
also specify `particle_type = photon`. We then move into discretization parameters, specifying that
linear Lagrange basis functions should be used and that the angular quadrature set should use
10 polar angles per half of the unit interval. Afterwards, we set `vacuum_boundaries = right` and
`reflective_boundaries = left` to apply our boundary conditions. Finally, we add volumetric sources to
block `0` (region 1) and block `3` (region 4). These sources are isotropic, and so we set `anisotropy` to
zero for both. We note that `volumetric_source_moments` is a double vector, where the first index corresponds
to the source block, and the second index is the source coefficients for a spherical harmonics expansion
in angle for all energy groups.

After spetting up our transport system, we proceed to add the material properties necessary for a transport
calculation. This is done with the [TransportMaterials block](AddTransportMaterialAction.md), which links
the added material properties to a transport system.

!listing /tutorials/1D_reed_sn/1D_reed.i
  block=TransportMaterials

We add a [ConstantTransportMaterial](ConstantTransportMaterial.md) for each region, which allows for the specification
of constant material properties in a spatial region. `transport_system = Neutron` links these properties to
the transport system we created earlier, and `block` specifies the subdomain that the materials should be
applied over. `group_total` is equivalent to $\Sigma_{t,g}$ (the group-wise total cross section), and
`group_scattering` is equivalent to $\Sigma_{s,g'\rightarrow g, l}$ (the group-wise scattering matrix).

The final components of a Gnat input file are the specification of the execution type (`Steady` vs `Transient` vs `Eigen`)
and the execution settings. The Reed problem is a fixed source problem, so we select `type = Steady` in the
`Executioner` block. We chose a PJFNK solver to avoid the need to form the full system Jacobian matrix which
results in substantial memory savings. The remaining settings improve the convergence rate of the PJFNK solver.
These execution settings can be found in the code block below.

!listing /tutorials/1D_reed_sn/1D_reed.i
  block=Executioner

In addition to specifying execution settings, we also add a vector post-processor to plot the solution
along the 1D mesh and specify that we want both Exodus and CSV output. This allows us to visualize our
results.

!listing /tutorials/1D_reed_sn/1D_reed.i
  start=VectorPostprocessors

## Execution and Post-Processing

To run the Reed benchmark problem,

```bash
gnat-opt -i 1D_reed.i
```

Congratulations, you've successfully run your first transport problem with Gnat! We can plot the output
from the vector post-processor (`1D_reed_out_scalar_flux_0001.csv`) in Python, which results in [reed_coarse].

!media media/1D_reed/coarse.png id=reed_coarse caption=Flux from the Reed benchmark with no mesh subdivisions.
  style=width:80%;margin-left:auto;margin-right:auto;halign:center

Our flux solution is quite coarse, indicating that added mesh refinement is required to resolve the spatial gradients.
The mesh can be refined by adding the following command line arguement when running the input file:

```bash
gnat-opt -i 1D_reed.i Mesh/uniform_refine=5
```

We can see that this took slightly longer to run then the coarse calculations (not by much as 1D transport problems are
quite fast); the results of which can be found in [reed_fine].

!media media/1D_reed/fine.png id=reed_fine caption=Flux from the Reed benchmark with 5 uniform subdivisions.
  style=width:80%;margin-left:auto;margin-right:auto;halign:center
