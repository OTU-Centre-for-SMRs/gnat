!config navigation breadcrumbs=False scrollspy=False

# HOME style=visibility:hidden;

# Gnat class=center light style=font-size:300%

# An open-source Moose application for the simulation of neutron activation
  class=center
  style=font-weight:200;font-size:200%

Gnat (+G+eneric +n+eutron +a+ctivation +t+oolkit) is a Moose-based finite
element application to study the behaviour of materials which have been
activated in a neutron field, serving as a source-term solver for Caribou. Gnat
provides a +S@n@+ neutron transport solver geared towards source-driven and
transient radiation fields encountered in neutron activation and shielding
scenarios. Gnat also provides a mass transport model for the activation of
species in both liquids and gases and a depletion solver for activated
materials. Moose provides several advantages for Gnat: scalability to machines
from single-digit core counts to tens of thousands of cores, continuous
open-source framework improvements, semi-automatic software quality assurance
(to NQA-1), and additional physics modules (such as finite volume and finite
  element fluid mechanics).

!row!
!col! small=12 medium=6 large=6 icon=device_hub

### Design and Theory class=center style=font-weight:200;

Familiarize yourself with the [design](syntax/NeutronActivationStudy/index.md)
of Gnat and the governing equations implemented:

- [Governing Equations](about/equations.md)
- [Equation Stabilization](about/stabilization.md)
- [Neutron Transport Angular Discretization](about/nte_angular_approach.md)
- [Numerical Implementation in Moose](about/implementation.md)

!col-end!

!col! small=12 medium=6 large=6 icon=school

### Installation and Examples class=center style=font-weight:200;

Get started with Gnat by installing the application and working through the
tutorials and example problems:

- [Installation](getting_started/installation.md)
- Tutorial 1: A Simple Neutron Transport Simulation
- Tutorial 2: The Duct Leg Problem
- Tutorial 3: Isotope Mass Transport
- Tutorial 4: Conjugate Neutron Transport and Activation
- [Other Examples](getting_started/other_examples.md)

!col-end!
!row-end!
