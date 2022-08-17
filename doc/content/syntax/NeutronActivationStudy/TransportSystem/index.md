# TransportSystem System

## Overview

The TransportSystem is responsible for initializing Gnat's built-in neutron
transport solver. Gnat discretizes the multi-group neutron transport equation (NTE)
with the discrete ordinates method in angle and applies the spherical harmonics
method to treat the scattering operator. This results in a large number of
nonlinear unknowns in the form of angular flux ordinates, and a number of
auxiliary variables in the form of moments of the angular neutron flux. In
addition to the sheer number of variables generated due to the angular
discretization, the multi-group formalism also yields coupled non-linear
variables in many energy groups. To combat the number of MooseObjects required to
solve this form of the NTE, the TransportSystem block initializes a single
[NeutronTransportAction](source/actions/NeutronTransportAction.md) which is
responsible for performing the following:

- Adding nonlinear variables for the angular flux ordinates for each energy group.
- Adding auxiliary variables for the moments of the angular flux for each energy group.
- Determining the spatial discretization and stabilization scheme from user input.
- Adding kernels (normal, discontinuous, and Dirac) required by the chosen scheme for the angular flux ordinates.
- Coupling the different group angular flux ordinates together through the scattering operator.
- Adding boundary conditions specified by the user for the angular flux ordinates.
- Adding initial conditions for the flux ordinates.
- Adding auxiliary kernels required to compute the angular flux moments.
- Adding custom output objects to minimize exodus file clutter.

## Constructed MooseObjects

[Table 1](#transport-objects) is a list of MooseObjects constructed by the
TransportSystem. Objects which contain SAAF correspond to the self-adjoint angular flux method while objects that contain DFEM correspond to the upwinding discontinuous Galerkin method.

!table id=transport-objects caption=MooseObjects constructed by the TransportSystem.
| Name | Type | Functionality | Condition(s) | Number Added |
| - | - | - | - | - |
| ADSAAFTimeDerivative | ADKernel | Applies the SAAF-stabilized time derivative term in the NTE to a single flux ordinate. | scheme = saaf_cfem,\\ execution_type = transient | $N\times G$ |
| ADSAAFStreaming | ADKernel | Applies the SAAF-stabilized streaming term (with a void treatment) in the NTE to a single flux ordinate. | scheme = saaf_cfem | $N\times G$ |
| ADSNRemoval | ADKernel | Applies the collision term in the NTE to a single flux ordinate. | N/A | $N\times G$ |
| ADSAAFMaterialSource | ADKernel | Applies the SAAF-stabilized external source term in the NTE to a single flux ordinate. Requires external source moments provided by the material system. | scheme = saaf_cfem | $N\times G$ |
| ADSAAFExternalScattering | ADKernel | Applies the SAAF-stabilized out of group scattering source term in the NTE to a single flux ordinate. Requires flux moments from the previous scattering source iteration. | scheme = saaf_cfem,\\ num_groups > 1,\\ debug_disable_source_iteration = false,\\ debug_disable_scattering = false | $N\times G$ |
| ADSAAFInternalScattering | ADKernel | Applies the SAAF-stabilized within-group scattering source term in the NTE to a single flux ordinate. Requires flux moments from the previous scattering source iteration. | scheme = saaf_cfem,\\ debug_disable_source_iteration = false,\\ debug_disable_scattering = false | $N\times G$ |
| ADSAAFScattering | ADKernel | Applies the SAAF-stabilized complete scattering source term in the NTE to a single flux ordinate. Requires all angular flux ordinates. | scheme = saaf_cfem,\\ debug_disable_source_iteration = true,\\ debug_disable_scattering = false | $N\times G$ |
| ADDFEMTimeDerivative | ADKernel | Applies the first order time derivative term in the NTE to a single flux ordinate. | scheme = upwinding_dfem,\\ execution_type = transient | $N\times G$ |
| ADDFEMStreaming | ADKernel | Applies the first order streaming term in the NTE to a single flux ordinate. | N/A | $N\times G$ |
| ADDFEMMaterialSource | ADKernel | Applies the first order external source term in the NTE to a single flux ordinate. Requires external source moments provided by the material system. | scheme = upwinding_dfem | $N\times G$ |
| ADDFEMExternalScattering | ADKernel | Applies the first order out of group scattering source term in the NTE to a single flux ordinate. Requires flux moments from the previous scattering source iteration. | scheme = upwinding_dfem,\\ num_groups > 1,\\ debug_disable_source_iteration = false,\\ debug_disable_scattering = false | $N\times G$ |
| ADDFEMInternalScattering | ADKernel | Applies the first order within-group scattering source term in the NTE to a single flux ordinate. Requires flux moments from the previous scattering source iteration. | scheme = upwinding_dfem,\\ debug_disable_source_iteration = false,\\ debug_disable_scattering = false | $N\times G$ |
| ADDFEMScattering | ADKernel | Applies the first order complete scattering source term in the NTE to a single flux ordinate. Requires all angular flux ordinates. | scheme = upwinding_dfem,\\ debug_disable_source_iteration = true,\\ debug_disable_scattering = false | $N\times G$ |
| SAAFIsoPointSource | DiracKernel | Applies the SAAF-stabilized external source term in the NTE to a single flux ordinate as an isotropic point source. Requires a list of `l = 0` source moments, a list of point source locations, and the point source group emission spectra. | scheme = saaf_cfem,\\ $n_{src}$ > 0 | $N\times n_{src}$ |
| DFEMIsoPointSource | DiracKernel | Applies the first order external source term in the NTE to a single flux ordinate as an isotropic point source. Requires a list of `l = 0` source moments, a list of point source locations, and the point source group emission spectra. | scheme = upwinding_dfem,\\ $n_{src}$ > 0 | $N\times n_{src}$ |
| ADDFEMUpwinding | ADDGKernel | Applies the first-order upwind stabilization term on element faces for the discontinuous Galerkin finite element method. | scheme = upwinding_dfem | $N\times G$ |
| NeutronFluxMoment | AuxKernel | Computes the moment of the angular neutron flux of degree `l`, order `m` for the neutron group `g`. Requires all flux ordinates for the current group. | N/A | $\begin{cases} G\times (L_{req} + 1)& 1D \\G\times 0.5(L_{req} + 1)(L_{req} + 2) & 2D\\G\times (L_{req} + 1)^{2} & 3D \end{cases}$ |
| ADSNVacuumBC | ADIntegratedBC | Applies the vacuum boundary condition in the NTE to a single flux ordinate. | $n_{v}$ > 0 | $n_{v}\times N\times G$ |
| ADSNMatSourceBC | ADIntegratedBC | Applies the surface source boundary condition in the NTE to a single flux ordinate. Requires surface source moments provided by the material system. | $n_{s}$ > 0 | $n_{s}\times N\times G$ |
| ADSNReflectiveBC | ADIntegratedBC | Applies the reflective boundary condition in the NTE to a single flux ordinate. | $n_{r}$ > 0 | $n_{r}\times N\times G$ |
| ConstantIC | InitialCondition | Initializes all angular flux ordinates with a constant value provided per energy group. | execution_type = transient,\\ ic_type = constant | $N\times G$ |
| Exodus | Output | Outputs the results of the neutron transport simulation to an exodus file with the angular neutron fluxes pruned if the user requests. | N/A | 1 |

$N$ is the number of flux ordinates.\\
$G$ is the number of energy groups.\\
$n_{src}$ is the minimum of the number of elements provided in `point_source_locations`,  `point_source_intensities` and `point_source_groups`.\\
$L_{req}$ is the maximum moment degree that Gnat should evaluate, as requested by the user, and corresponds to `max_anisotropy`.\\
$n_{v}$, $n_{s}$ and $n_{r}$ are the number of boundaries with vacuum, surface source and reflective boundary conditions (respectively).

## Required Material Properties id=required-materials

The MooseObjects constructed by the TransportSystem require a large number of
material properties, some of which are only required for specific circumstances.
A list of properties can be found below in [Table 2](#transport-materials):

!table id=transport-materials caption=Material properties required by the TransportSystem MooseObjects.
| Property Name | Type | Description | Consuming MooseObjects | Required for | Number Required |
| - | - | - | - | - | - |
| v_g | std::vector<Real> | The average neutron speed in all groups. | ADSAAFTimeDerivative\\ ADDFEMTimeDerivative  | Transient simulations. | $G$ |
| removal_xs_g | std::vector<Real> | The neutron removal cross-section in all groups. | ADSNRemoval\\ ADSAAFTimeDerivative \\ ADSAAFStreaming \\ ADSAAFMaterialSource \\ ADSAAFExternalScattering \\ ADSAAFInternalScattering \\ ADSAAFScattering \\ SAAFIsoPointSource | All. | $G$ |
| surface_source | std::vector<Real> | Surface sources for all flux ordinates and groups. | ADSNMatSourceBC | Surface source boundary conditions. | $N\times G$ |
| scattering_matrix | std::vector<Real> | The scattering matrix where each entry is composed of scattering cross-section moments. | ADSAAFExternalScattering \\ ADSAAFInternalScattering \\ ADSAAFScattering \\ ADDFEMExternalScattering \\ ADDFEMInternalScattering \\ ADDFEMScattering | Simulations with scattering. | $L_{scat}\times G\times G$ |
| medium_anisotropy | unsigned int | The scattering anisotropy (degree) of the material. | ADSAAFExternalScattering \\ ADSAAFInternalScattering \\ ADSAAFScattering \\ ADDFEMExternalScattering \\ ADDFEMInternalScattering \\ ADDFEMScattering | Simulations with scattering. | 1 |
| source_moments | std::vector<Real> | The external volumetric source moments for all groups. | ADSAAFMaterialSource\\ ADDFEMMaterialSource | Simulations with volumetric sources. | $\begin{cases} G\times (L_{src} + 1)& 1D \\G\times 0.5(L_{src} + 1)(L_{src} + 2) & 2D\\G\times (L_{src} + 1)^{2} & 3D \end{cases}$ |
| medium_source_anisotropy | unsigned int | The external volumetric source anisotropy (degree) of the material. | ADSAAFMaterialSource\\ ADDFEMMaterialSource | Simulations with volumetric sources. | 1 |
| saaf_eta | Real | The SAAF stabilization cutoff. | ADSAAFTimeDerivative \\ ADSAAFStreaming \\ ADSAAFMaterialSource \\ ADSAAFExternalScattering \\ ADSAAFInternalScattering \\ ADSAAFScattering \\ SAAFIsoPointSource | SAAF simulations. | 1 |
| saaf_c | Real | The SAAF stabilization parameter. | ADSAAFTimeDerivative \\ ADSAAFStreaming \\ ADSAAFMaterialSource \\ ADSAAFExternalScattering \\ ADSAAFInternalScattering \\ ADSAAFScattering \\ SAAFIsoPointSource | SAAF simulations. | 1 |

$L_{scat}$ is the maximum anisotropy (degree) of the scattering material
corresponding to `medium_anisotropy`.\\
$L_{src}$ is the maximum anisotropy (degree) of the external source
corresponding to `medium_source_anisotropy`.

Gnat provides several material by default. An up-to-date list of these materials
can be found on the [Gnat sources page](source/index.md).

## Example Input File Syntax

There are two ways to employ the TransportSystem block. The first is to use it
in conjunction with the Common Parameters block:

```language=moose
[NeutronActivationStudy]
  execution_type = steady
  num_groups = 1
  max_anisotropy = 0

  debug_verbosity = level0

  [TransportSystem]
    scheme = saaf_cfem
    output_angular_fluxes = true

    order = FIRST
    family = LAGRANGE

    n_azimuthal = 2
    n_polar = 2

    vacuum_boundaries = 'left right top bottom'

    point_source_locations = '5.0 5.0 0.0'
    point_source_intensities = '1000.0'
    point_source_groups = '1'
  []
[]
```

[/examples/2D_scattering_cgfem/test.i]

The other option is to avoid the Common Parameters block. This can be useful if
only the neutron transport capabilities of Gnat are required:

```language=moose
[NeutronActivationStudy/TransportSystem]
  scheme = saaf_cfem
  execution_type = steady
  num_groups = 1
  max_anisotropy = 0

  debug_verbosity = level0
  output_angular_fluxes = true

  order = FIRST
  family = LAGRANGE

  n_azimuthal = 2
  n_polar = 2

  vacuum_boundaries = 'left right top bottom'

  point_source_locations = '5.0 5.0 0.0'
  point_source_intensities = '1000.0'
  point_source_groups = '1'
[]
```

Both of these approaches yield the same results. However, the second approach will not
auto-fill required parameters for transport materials.

The example provided above initializes a steady-state (`execution_type = steady`)
single group (`num_groups = 1`) transport simulation with the SAAF method
(`scheme = saaf_cfem`). Vacuum boundaries are specified for all sides. An
isotropic point source is specified at `(5.0, 5.0, 0.0)` with an $l = 0$ moment
of `1000` emitting into energy group `1`.

!alert note
The choice of a transport scheme must be paired with an appropriate set of
finite element basis functions. The SAAF method must use continuous basis
functions and the upwinding method must use discontinuous basis functions.

The following example is the same sceneario as described above, but with the
upwinding scheme and discontinuous basis functions:

```language=moose
[NeutronActivationStudy]
  execution_type = steady
  num_groups = 1
  max_anisotropy = 0

  debug_verbosity = level0

  [TransportSystem]
    scheme = upwinding_dfem
    output_angular_fluxes = true

    order = FIRST
    family = MONOMIAL

    n_azimuthal = 2
    n_polar = 2

    vacuum_boundaries = 'left right top bottom'

    point_source_locations = '5.0 5.0 0.0'
    point_source_intensities = '1000.0'
    point_source_groups = '1'
  []
[]
```

[/examples/2D_scattering_dgfem/test.i]

A final example is the case of multi-group neutron transport with the TransportSystem action:

```language=moose
[NeutronActivationStudy]
  execution_type = steady
  num_groups = 2
  max_anisotropy = 0

  [TransportSystem]
    scheme = saaf_cfem
    output_angular_fluxes = false

    order = FIRST
    family = LAGRANGE

    n_azimuthal = 2
    n_polar = 2

    vacuum_boundaries = 'left right top bottom'

    point_source_locations = '5.0 5.0 0.0'
    point_source_intensities = '1000.0'
    point_source_groups = '1'
  []
[]
```

[/examples/2D_scattering_2_group_cgfem/test.i]

Two neutron groups are used (hence `num_groups = 2`). It's possible to have an
isotropic point source which emits into multiple groups:

```language=moose
point_source_locations = '5.0 5.0 0.0; 5.0 5.0 0.0'
point_source_intensities = '1000.0 1000.0'
point_source_groups = '1 2'
```

The above example source emits equally in both groups.

!alert warning
The isotropic point source syntax is likely to change in the future as the
current implementation is quite cumbersome.

!syntax list /NeutronActivationStudy/TransportSystem objects=True actions=False subsystems=False

!syntax list /NeutronActivationStudy/TransportSystem objects=False actions=False subsystems=True

!syntax list /NeutronActivationStudy/TransportSystem objects=False actions=True subsystems=False
