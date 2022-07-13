# TODO List for Gnat

## Neutron Transport Solver MVP

- Add support for 1D and 2D cartesian geometry.
  - Modification of the quadrature set class and SPH functions.
  - Reduced order quadrature sets and SPH moments.

- ConstantNeutronicsMaterial:
  - Isotropic scattering and source, other properties are constant.

- Add the NeutronTransportAction class. Sets up the transport problem.

## Mass Transport and Activation Solver MVP

### Fluids

- Kernels
  - ADIsotopeTimeDerivative
  - ADIsotopeAdvection : ConservativeAdvection
  - ADIsotopeDiffusion
  - ADIsotopeFormation (takes in microscopic cross-sections)
  - ADIsotopeDecay

- DGKernels
  - ADDGIsotopeUpwinding

- BCs
  - ADIsotopeOutflowBC
  - ADIsotopeInflowBC

### Solids

- ODEs
  - Bateman equations for activation and removal. Solve with CRAM?

### Both

- Materials
  - IsotopeMaterial : microscopic cross-sections + density (no coupled density variable), or just microscopic cross-sections (coupled density variable). Decay constant.

- Actions
  - NeutronActivationAction : solids and fluids handled separately.

## Cross-section generation?
- Figure this out I guess.
