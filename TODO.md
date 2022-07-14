# TODO List for Gnat

## Neutron Transport Solver MVP

- Add support for 1D and 2D cartesian geometry.
  - Modification of the quadrature set class and SPH functions.
    - Major axis of 2D problems is the x-axis. Only have to evaluate moments where
    m > 0 and quadrature set vectors with _direction(2) > 0 since the problem is
    symmetrical about the x-y plane. Need to multiply expanded sources
    (external and scattering) by 2 to preserve flux symmetry.
    - Major axis of 1D problems is naturally the x-axis. Only need Legendre
    moments of the flux and external source. Need to multiply expanded sources
    (external and scattering) by 2\pi to preserve flux symmetry.

- ConstantNeutronicsMaterial:
  - Isotropic scattering and source, other properties are constant.

- Add the NeutronTransportAction class. Sets up the transport problem.
  - Specify quadrature orders as direction cosines per quadrant. Makes quadrature
  sets consistent between 1D, 2D and 3D problems.

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
