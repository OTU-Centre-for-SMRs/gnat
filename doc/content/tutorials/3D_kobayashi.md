# Modelling the Kobayashi Benchmark Problems

In this tutorial, you will learn how to model more complicated 3D fixed source problems in [!cite](kobayashi) with Gnat.
To access this tutorial,

```bash
cd gnat/tutorials/3D_kobayashi
```

## Common Material Properties

All of the Kobayashi problems use the same material properties. These consist of an isotropic neutron
source, a near-void material (to generate streaming paths), and a shield. The specifications for the
materials can be found in [kob_mats].

!table id=kob_mats caption=Material properties for the Kobayashi benchmarks.
| Material | $\Sigma_{t}$ ($\text{cm}^{-1}$) | $\Sigma_{s}$ ($\text{cm}^{-1}$) | $q$ ($\text{cm}^{-3}s^{-1}$) |
| :- | - | - | - |
| Source | $1\times10^{-1}$ | $5\times10^{-2}$ | 1.0 |
| Near-Void | $1\times10^{-4}$ | $5\times10^{-5}$ | 0.0 |
| Shield | $1\times10^{-1}$ | $5\times10^{-2}$ | 0.0 |

## Volume Source in a Cubic Shield

The first Kobayashi benchmark problem consists of a $10\text{ cm }\times 10\text{ cm }\times 10\text{ cm}$
cubic volume source with one corner centered on the origin. The source is surrounded on three sides by a
$40\text{ cm}$ thick near-source reigon, which is further encased in a $50\text{ cm}$ thick shield.
The boundary planes with $x = 0$, $y = 0$, and $z = 0$ use reflective boundary conditions. The boundary
planes with $x = 100$, $y = 100$, and $z = 100$ are vaccum boundaries. The geometry of the problem can be
found in [kob_1_geom].

!media media/3D_kobayashi/kob_1_geom.png id=kob_1_geom caption=Geometry of the first Kobayashi problem. Green: source. White: near-void. Red: shield.
  style=width:50%;margin-left:auto;margin-right:auto;halign:center

The Gnat input file to model this problem starts with the mesh, which is generated with a
[CartesianMeshGenerator](https://mooseframework.inl.gov/source/meshgenerators/CartesianMeshGenerator.html).
We specify that the generated mesh should be 3D (`dim = 3`) and that we have three distinct subdivisions with
`dx`, `dy`, and `dz`. We then set the number of elements we want for each subdivision with `ix`, `iy`, and
`iz`. In this model, we start with one element per $10\text{ cm}$ interval. We then add block ids for material
and volume source assignment using `subdomain_id`. This vector must be the same length as the product of the
lengths of `dx`, `dy`, and `dz`. The generated mesh is uniformely refined once with `uniform_refine = 1`.
This mesh generator automatically adds side sets for the six faces of the cube named `left`, `right`, `top`,
`bottom`, `front`, and `back` which we can use for boundary conditions.

!listing /tutorials/3D_kobayashi/kobayashi_1.i
  block=Mesh

After specifying the mesh, we move on to setting up the transport system. The Kobayashi problems are mono-energetic,
so we specify that we only want a single energy group. Similarily, the scattering cross sections are isotropic
and so we don't need to evaluate higher order angular moments of the flux. We specify that we want to use the
SAAF-CFEM-S$_{\text{N}}$ method with `scheme = saaf_cfem`, and that `particle_type = neutron` (which is arbitrary,
`particle_type = photon` would also suffice). We specify that linear Lagrange basis functions should be used, and that
a Gauss-Chebyshev angular quadrature set with five polar angles and five azimuthal angles per octant of the unit sphere
should be used. Finally, we apply reflective boundary conditions along the lines of symmetry, apply vacuum boundary
conditions to the remaining faces, and specify that the isotropic volume source with unit intensity should be applied
to block 2.

!listing /tutorials/3D_kobayashi/kobayashi_1.i
  block=TransportSystems

Next, the material properties of the domain must be specified with `TransportMaterials`:

!listing /tutorials/3D_kobayashi/kobayashi_1.i
  block=TransportMaterials

This declares three a constant materials which provides enough parameters for a single
neutron group. `Shield` implements the shield material and souce material from [kob_mats],
while `Duct` implements the near-void material. We note that the main differences between
`Shield` and `Duct` are the values of `group_total`, `group_scattering`, and `block`. Both materials
set `transport_system = 'Neutron'` to link to the single transport system we set up, and specify
`anisotropy = 0` (as these provided scattering cross sections are isotropic).

The remainder of the input file describes how the problem should be solved, the post-processors
that should be added, and the outputs requested.

!listing /tutorials/3D_kobayashi/kobayashi_1.i
  block=Executioner

The `Executioner` block declares the solver parameters and solve type. In this case,
`type = Steady` and the Preconditioned Jacobian-Free Newton Krylov
(`solve_type = PJFNK`) is used. It is recommended that the user employ the PJFNK
solver for continuous finite element methods when solving radiation transport
problems with the SAAF-CFEM-S$_{\text{N}}$ scheme. Next, we add a collection of
`VectorPostprocessors` to plot the solution along the lines required to match the results
reported by [!cite](kobayashi). These include: a line from the center of the source
along the y direction (`Line_A`); a line along the diagonal from the center of the source to
the corner of the domain (`Line_B`); and a line along the discontinuity between the shield
and the near-void region (`Line_C`).

!listing /tutorials/3D_kobayashi/kobayashi_1.i
  block=VectorPostprocessors

As we are adding post-processors, we need to specify that we want CSV output in addition
to Exodus output.

!listing /tutorials/3D_kobayashi/kobayashi_1.i
  block=Outputs

The simulation can be executed with the following shell command:

```language=bash
mpiexec -np 4 gnat-opt -i ./kobayashi_1.i --n-threads=2
```

which will run the problem with four MPI ranks and two OpenMP threads per rank. The scalar flux
calculated by Gnat can be found in [kob_1_flux].

!media media/3D_kobayashi/kob_1_flux.png id=kob_1_flux caption=Scalar flux for the first Kobayashi benchmark.
  style=width:50%;margin-left:auto;margin-right:auto;halign:center

Plotting the post-processed results from the three [LineValueSampler](https://mooseframework.inl.gov/source/vectorpostprocessors/LineValueSampler.html)s
results in [kob_1a], [kob_1b], and [kob_1c].

!row! style=display:inline-flex;
!col!

!media media/3D_kobayashi/kob_1_A_err.png id=kob_1a caption=Scalar flux plotted along line A for the first Kobayashi benchmark.
  style=width:100%;margin-left:auto;margin-right:auto;halign:center

!col-end!

!col!

!media media/3D_kobayashi/kob_1_B_err.png id=kob_1b caption=Scalar flux plotted along line B for the first Kobayashi benchmark.
  style=width:100%;margin-left:auto;margin-right:auto;halign:center

!col-end!

!col!

!media media/3D_kobayashi/kob_1_C_err.png id=kob_1c caption=Scalar flux plotted along line C for the first Kobayashi benchmark.
  style=width:100%;margin-left:auto;margin-right:auto;halign:center

!col-end!

!row-end!

We can clearly see the impact of ray effects from the limited number of samples in the angular
quadrature set, alongside the imapct of the coarse spatial discretization which manifest in
oscillations in the scalar flux. These will be come more pronounced in the following sections,
where the straight duct and dog-legged duct problems are modelled.

## Straight Duct

The second Kobayashi benchmark problem consists of a $10\text{ cm }\times 10\text{ cm }\times 10\text{ cm}$
cubic volume source with one corner centered on the origin. It is surrounded by a shield that is
$90\text{ cm}$ thick along the positive y direction and $50\text{ cm}$ thick along the x and z direction.
The $90\text{ cm}$ long component of the shield is penetrated by a duct with a square cross section,
yielding a long streaming path. The geometry of this problem can be found in [kob_2_geom]. Similar
to the first benchmark problem, three sides use reflective boundary conditions while the remaining
three sides use vacuum boundaries.

!media media/3D_kobayashi/kob_2_geom.png id=kob_2_geom caption=Geometry of the second Kobayashi problem. Green: source. White: near-void. Red: shield.
  style=width:50%;margin-left:auto;margin-right:auto;halign:center

There are two differences between the previously modelled cubic shield problem and the straight
duct problem. The first is the definition of the mesh where `dx`, `dy`, and `dz` are modified to
change the configuration of the geometry. `ix`, `iy`, and `iz` are then changed to maintain
$10\text{ cm}$ element intervals. Finally, the blocks are remapped to capture the material changes.
The changes to the mesh can be found below.

!listing /tutorials/3D_kobayashi/kobayashi_2.i
  block=Mesh

The second change to the input file is the [LineValueSampler](https://mooseframework.inl.gov/source/vectorpostprocessors/LineValueSampler.html)s
we use to post-process the solution. `Line_A` is modified such that it plots down the length
of the duct and `Line_B` is modified such that it plots from the centerline of the duct
into the shield at the edge of the domain. The new `VectorPostprocessor`s can be found below.

!listing /tutorials/3D_kobayashi/kobayashi_2.i
  block=VectorPostprocessors

The second Kobayashi problem can be run with the following:

```language=bash
mpiexec -np 4 gnat-opt -i ./kobayashi_2.i --n-threads=2
```

The scalar flux calculated by Gnat can be found in [kob_2_flux].

!media media/3D_kobayashi/kob_2_flux.png id=kob_2_flux caption=Scalar flux for the second Kobayashi benchmark.
  style=width:50%;margin-left:auto;margin-right:auto;halign:center

We can really see the ray effects in the solution which are caused by the low order quadrature
set. Increasing the number of azimuthal angles per octant (`n_azimuthal`) would aid in mitigating
these ray effects by adding more quadrature directions along the duct, though it would result in
a large runtime penalty. The error between the Gnat results and the benchmark along `Line_A` and
`Line_B` are plotted in [kob_2a] and [kob_2b]. Oscillations caused by the coarse mesh can be seen
in [kob_2a], and large undershoots in the solution compared to the reference (due to ray effects)
can be seen in [kob_2b].

!row! style=display:inline-flex;
!col!

!media media/3D_kobayashi/kob_2_A_err.png id=kob_2a caption=Scalar flux plotted along line A for the second Kobayashi benchmark.
  style=width:100%;margin-left:auto;margin-right:auto;halign:center

!col-end!

!col!

!media media/3D_kobayashi/kob_2_B_err.png id=kob_2b caption=Scalar flux plotted along line B for the second Kobayashi benchmark.
  style=width:100%;margin-left:auto;margin-right:auto;halign:center

!col-end!

!row-end!

## Dog-Legged Duct

The last Kobayashi problem is similar to the second problem, where the only differece is that
the duct makes three $90^{\circ}$ turns before exiting the shielding medium. The geometry of this
problem can be found in [kob_3_geom_1] and [kob_3_geom_2].

!row! style=display:inline-flex;
!col!

!media media/3D_kobayashi/kob_3_geom_1.png id=kob_3_geom_1 caption=Full geometry of the third Kobayashi problem. Green: source. White: near-void. Red: shield.
  style=width:100%;margin-left:auto;margin-right:auto;halign:center

!col-end!

!col!

!media media/3D_kobayashi/kob_3_geom_2.png id=kob_3_geom_2 caption=Dog-legged duct in the third Kobayashi problem.
  style=width:100%;margin-left:auto;margin-right:auto;halign:center

!col-end!

!row-end!

Similar to the straight duct, the only changes we need to make to the input file are to the mesh and
the post-processors we add to compare the solution to the reference solution. The `[Mesh]` block for
this problem is substantially more complicated then the first and second problem due to the added
heterogeneity; extra slices in `dx`, `dy`, and `dz` are required to model these material changes alongside
additional block mappings in `subdomain_id`. The full mesh block can be found below.

!listing /tutorials/3D_kobayashi/kobayashi_3.i
  block=Mesh

When it comes to the changes in post-processors, `Line_A` plots along the
center of the duct (before the first bend). `Line_B` plots along the interface between the duct and the
shield at the first bend. `Line_C` plots across the duct after the third bend. The modifiications to the
`LineValueSampler` post-processors can be found below.

!listing /tutorials/3D_kobayashi/kobayashi_3.i
  block=VectorPostprocessors

The final Kobayashi problem can be run with the following:

```language=bash
mpiexec -np 4 gnat-opt -i ./kobayashi_3.i --n-threads=2
```

The scalar flux over the whole problem can be found in [kob_3_flux_1], and the scalar flux along the duct
can be found in [kob_3_flux_2].

!row! style=display:inline-flex;
!col!

!media media/3D_kobayashi/kob_3_flux_1.png id=kob_3_flux_1 caption=Scalar flux in the third Kobayashi benchmark.
  style=width:100%;margin-left:auto;margin-right:auto;halign:center

!col-end!

!col!

!media media/3D_kobayashi/kob_3_flux_2.png id=kob_3_flux_2 caption=Scalar flux isolated to the dog-legged duct.
  style=width:100%;margin-left:auto;margin-right:auto;halign:center

!col-end!

!row-end!

The results along each line compared to the reference scalar flux can be found below.

!row! style=display:inline-flex;
!col!

!media media/3D_kobayashi/kob_3_A_err.png id=kob_3a caption=Scalar flux plotted along line A for the first Kobayashi benchmark.
  style=width:100%;margin-left:auto;margin-right:auto;halign:center

!col-end!

!col!

!media media/3D_kobayashi/kob_3_B_err.png id=kob_3b caption=Scalar flux plotted along line B for the first Kobayashi benchmark.
  style=width:100%;margin-left:auto;margin-right:auto;halign:center

!col-end!

!col!

!media media/3D_kobayashi/kob_3_C_err.png id=kob_3c caption=Scalar flux plotted along line C for the first Kobayashi benchmark.
  style=width:100%;margin-left:auto;margin-right:auto;halign:center

!col-end!

!row-end!
