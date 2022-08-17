# NeutronActivationStudy System

## Overview

The neutron activation system is responsible for initializing and setting up
steady-state or transient neutron activation problems in Gnat. The system is
composed of three components:

- The common parameters list.
- The TransportSystem block.
- The ActivationSystem block.

Under the hood, the different blocks translate to MOOSE actions which initialize
the many MOOSE objects required to set up neutron transport and activation
studies.

### The Common Parameters List

The common parameters list provides parameters necessary for both neutron
transport simulations and neutron activation studies. This avoids the need to
declare the same common parameters multiple times across a single Gnat input file.
The common parameters list aliases the construction of a
[CommonGnatAction](source/actions/CommonGnatAction.md). An example of how to use
the common parameters list can be found below:

```language=moose
[NeutronActivationStudy]
  # The common parameters list:
  num_groups = 1 # The number of spectral energy groups for neutron transport simulations.
  execution_type = steady # The execution type for the problem. Either steady or transient.

  # Other blocks go here afterwards...
  [TransportSystem]
    # ...
  []

  [ActivationSystem]
    # ...
  []
[]
```

The parameters `num_groups` and `execution_type` are mandatory if the user
plans on using Gnat's built-in neutron transport capabilities. A partial list of
common parameters can be found in [Table 1](#common-params), and the full list
in [CommonGnatAction](source/actions/CommonGnatAction.md).

!table id=common-params caption=Common parameters block parameters.
| Name | Parameter Group | Value | Default | Description |
| - | - | - | - | - |
| execution_type | Required | MooseEnum | N/A | The execution type of the problem. Options are `steady` and `transient`. |
| num_groups | Required | unsigned int | N/A | The number of spectral energy groups. The provided number of groups must be greater than 0. |
| block | Optional | std::vector<SubdomainName> | N/A | The list of blocks (ids or names) that this variable will be applied. |
| debug_verbosity | Optional | MooseEnum | level1 | How verbose the debug output of the system should be. Options are `level0` (fully verbose) and `level1` (minimal). |
| flux_moment_names | Optional | std::string | flux_moment | Variable names for the moments of the angular flux. The output format for the group flux moments will be of the form `flux_moment_names`_g_l_m. |

### The TransportSystem Block

The TransportSystem is responsible for initializing the many MooseObjects
required to perform neutron transport simulations. While it is possible to
declare all of the required objects manually in a standard Moose input file the
number of unknowns in all transport schemes grows rapidly as the complexity of
the problem increases, rendering it impractical to do so. The TransportSystem
block aliases the construction of a
[NeutronTransportAction](source/actions/NeutronTransportAction.md). An example
of how to use the TransportSystem block can be found below:

```language=moose
[NeutronActivationStudy]
  # Common parameters go here...

  # The transport system block:
  [TransportSystem]
    scheme = saaf_cfem # The scheme to be employed for transport calculations.
    output_angular_fluxes = true # Whether or not to output the angular flux ordinates.

    order = FIRST # The order of the finite element shape functions which represent the angular flux ordinates and the moments of the angular flux.
    family = LAGRANGE # The family of the finite element shape functions which represent the angular flux ordinates and the moments of the angular flux.

    n_azimuthal = 2 # The number of azimuthal quadrature points per quadrant in the angular quadrature.
    n_polar = 2 # The number of polar quadrature points per quadrant in the angular quadrature.

    vacuum_boundaries = 'left right top bottom' # The locations where vacuum boundaries should be applied.

    point_source_locations = '5.0 5.0 0.0' # Isotropic point source locations.
    point_source_intensities = '1000.0' # Isotropic point source emission intensities.
    point_source_groups = '1' # The groups the point sources emit into.
  []

  # Other blocks go here afterwards...
  [ActivationSystem]
    # ...
  []
[]
```

The parameters `scheme`, `order`, and `family` are required when using Gnat's
built-in transport capabilities. If the common parameter block is not being
employed, both `num_groups` and `execution_type` must be provided separately.
A partial list of parameters for the TransportSystem can be found in
[Table 2](#transport-params), and the full list in
[NeutronTransportAction](source/actions/NeutronTransportAction.md).

!table id=transport-params caption=TransportSystem block parameters.
| Name | Parameter Group | Value | Default | Description |
| - | - | - | - | - | - |
| execution_type | Required\* | MooseEnum | N/A | The execution type of the problem. Options are `steady` and `transient`. \*If the common parameters block is used this parameter is no longer required. |
| scheme | Required | MooseEnum | N/A | The transport scheme to apply. Options are `saaf_cfem` and `upwinding_dfem`. |
| num_groups | Required\* | unsigned int | N/A | The number of spectral energy groups. The provided number of groups must be greater than 0. *If the common parameters block is used this parameter is no longer required. |
| order | Required | MooseEnum | FIRST | The order of the finite element shape functions to use for the angular flux ordinates and moments. |
| family | Required | MooseEnum | N/A | The family of finite element shape functions to use for the angular flux ordinates and moments. |
| angular_flux_names | Variable | std::string | angular_flux | Variable names for the angular flux. The output format for the group angular fluxes will be of the form `angular_flux_names`_g_n. |
| flux_moment_names | Variable | std::string | flux_moment | Variable names for the moments of the angular flux. The output format for the group flux moments will be of the form `flux_moment_names`_g_l_m. If the common parameters block is used this parameter will be overwritten. |
| output_angular_fluxes | Variable | bool | false | Whether the angular flux ordinates should be written to the exodus file or not. |
| scaling | Variable | std::string | 1.0 | A multiplicative scaling factor to apply to the flux ordinates and moments. |
| block | Simulation | std::vector<SubdomainName> | N/A | The list of blocks (ids or names) that this variable will be applied. If the common parameters block is used this parameter will be overwritten. |
| max_anisotropy | Simulation | unsigned int | 0 | The maximum degree of anisotropy to evaluate. Defaults to 0 for isotropic scattering. If the common parameters block is used this parameter will be overwritten. |
| ic_type | Initial Condition | MooseEnum | constant | The type of initial condition to use (if multiple are provided). Options are `constant`, `function` and `file`. |
| constant_ic | Initial Condition | std::vector<double> | 0.0 | A constant initial condition for the angular fluxes. The flux initial conditions must be listed in order of decreasing group energy. |
| reflective_boundaries | Boundary Condition | std::vector<BoundaryName> | N/A | The boundaries to apply reflective boundary conditions. |
| source_boundaries | Boundary Condition | std::vector<BoundaryName> | N/A | The boundaries to apply surface source boundary conditions. |
| vacuum_boundaries | Boundary Condition | std::vector<BoundaryName> | N/A | The boundaries to apply vacuum boundary conditions. |
| major_axis | Quadrature | MooseEnum | x | Major axis of the angular quadrature. Allows the polar angular quadrature to align with a cartesian axis with minimal heterogeneity. This parameter is only applied in 3D cartesian problems. |
| n_polar | Quadrature | unsigned int | 3 | Number of Legendre polar quadrature points in a single octant of the unit sphere. |
| n_azimuthal | Quadrature | unsigned int | 3 | Number of Chebyshev azimuthal quadrature points in a single octant of the unit sphere. |
| debug_verbosity | Debugging | MooseEnum | level1 | How verbose the debug output of the system should be. Options are `level0` (fully verbose) and `level1` (minimal). If the common parameters block is used this parameter will be overwritten. |
| debug_disable_source_iteration | Debugging | bool | true | Debug option to disable source iteration. |
| debug_disable_scattering | Debugging | bool | false | Debug option to disable scattering. |

### The ActivationSystem Block

!alert construction title=Incomplete Feature
The neutron activation capabilities in Gnat are still under heavy development.
This piece of documentation will change quite a bit as the system is completed.

## Example Input File Syntax

!alert warning title=Incomplete Feature
The neutron activation capabilities in Gnat are still under heavy development.
Examples are provided for the neutron transport capabilities in Gnat exclusively.

The following example sets up the steady-state neutron transport capabilities of
Gnat with both the TransportSystem and Common Parameters blocks.

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

!alert note
The use of Gnat's neutron transport capabilities requires material objects which
inherit from the [EmptyNeutronicsMaterial](source/materials/EmptyNeutronicsMaterial.md)
class defined along the same subdomain(s) as the transport simulation. Please see
the [materials section](syntax/NeutronActivationStudy/TransportSystem/index.md#required-materials)
in the TransportSystem documentation for more information.

!syntax list /NeutronActivationStudy objects=false actions=True subsystems=false

!syntax list /NeutronActivationStudy objects=false actions=false subsystems=True
