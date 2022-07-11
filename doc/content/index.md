!config navigation breadcrumbs=False scrollspy=False

# HOME style=visibility:hidden;

# Gnat class=center light style=font-size:300%

# An open-source MOOSE application for the simulation of neutron activation
  class=center
  style=font-weight:200;font-size:200%

Gnat (+G+eneric +n+eutron +a+ctivation +t+oolkit) is a MOOSE-based finite element application
to study the behaviour of materials which have been activated in a neutron field,
serving as a source-term solver for Caribou. Gnat provides a +S@n@+ neutron
transport solver geared towards source-driven and transient radiation fields encountered
in neutron activation and shielding scenarios. Gnat provides cross-section generation
capabilities through NJOY, a aerosol transport model for the activation of species
in both liquids and gases, and a depletion solver for activated materials. MOOSE
provides several advantages for Gnat: scalability to machines from single-digit
core counts to tens of thousands of cores, continuous open-source framework
improvements, semi-automatic software quality assurance (to NQA-1), and
additional physics modules (such as finite volume and finite element fluids).
