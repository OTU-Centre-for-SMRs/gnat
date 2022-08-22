# Numerical Implementation in Moose

Moose is a Galerkin finite element framework which provides several plug-and-play
systems to assist in the implementation of partial differential equations. Gnat
takes advantage of the kernel, auxiliary kernel, Dirac kernel, discontinuous
Galerkin kernel, boundary condition, material and action systems. The different
systems used by Gnat are summarized below:

!table id=used_moose_systems caption=Description of Moose systems used by Gnat.
| System Name | Description |
| - | - |
| Kernel | Evaluates the non-linear residual contribution of a term in the system of coupled partial differential equations for a single non-linear variable. |
| Auxiliary Kernel (AuxKernel) | Evaluates the value of a variable auxiliary to the problem.  |
| Dirac Kernel | Evaluates the non-linear residual contribution of a point source term in the system of coupled partial differential equations for a single non-linear variable. |
| Discontinuous Galerkin Kernel (DGKernel) | Evaluates the element interior face terms which occur when using the discontinuous finite element method. |
| Boundary Condition (BC) | Evaluates the non-linear residual contribution of the implicit boundary conditions in the system of coupled partial differential equations for a single non-linear variable. |
| Material | Provides domain-dependant material properties to other Moose systems. |
| Action | A system which allows Moose users to specialize and automate components of how Moose sets up, solves, and post-processes the problem.  |

## Neutron Transport Implementation

Each residual contribution in a partial differential equation can be assigned to
a specific Moose system. The [SAAF weak form](stabilization.md#saaf) of the
neutron transport equation with the
[+S@n@+ disretization](nte_angular_approach.md#sn_disc) applied in angle and
[+P@n@+ treatment](nte_angular_approach.md#scattering) applied on the scattering
operator is:

!equation id=saaf_sn_wf
\underbrace{\Big( \psi_{j} + \tau_{g}\vec{\nabla}\phi_{j}\cdot\hat{\Omega}_{n},\, \frac{\partial}{\partial t}\frac{\Psi^{h}_{g,n}}{v_{g}} \Big)_{V}}_{\text{Kernel}}
+ \underbrace{\Big\langle \psi_{j},\, \hat{n}\cdot\hat{\Omega}_{n}\Psi^{h}_{g,n} \Big\rangle_{\Gamma}}_{\text{Boundary Condition}}
+ \underbrace{\Big( \vec{\nabla}\psi_{j}\cdot\hat{\Omega}_{n},\, \tau_{g}\hat{\Omega}_{n}\cdot\vec{\nabla}\Psi_{g,n}^{h} + (\tau_{g}\Sigma_{r,\, g} - 1)\Psi_{g,n}^{h} \Big)_{V}}_{\text{Kernel}}\\
+ \underbrace{\Big( \psi_{j},\, \Sigma_{r,\,g}\Psi^{h}_{g,n} \Big)_{V}}_{\text{Kernel}}
- \underbrace{\Big( \psi_{j} + \tau_{g}\vec{\nabla}\phi_{j}\cdot\hat{\Omega}_{n},\, \sum_{g' = 1}^{G}\Sigma_{s,\, g'\rightarrow g}\sum_{l = 0}^{L}\frac{2l + 1}{4\pi} f_{g'\rightarrow g,\, l}\sum_{m = -l}^{l}Y_{l,m}(\hat{\Omega}_{n})\Phi_{g',l,m} \Big)_{V}}_{\text{Kernel}}
- \underbrace{\Big( \psi_{j} + \tau_{g}\vec{\nabla}\phi_{j}\cdot\hat{\Omega}_{n},\, S_{g,n} \Big)_{V}}_{\text{Kernel or Dirac Kernel}} = 0

Similarly for the [DGFEM upwinding weak form](stabilization.md#upwind) of the
neutron transport equation with the
[+S@n@+ disretization](nte_angular_approach.md#sn_disc) applied in angle and
[+P@n@+ treatment](nte_angular_approach.md#scattering) applied on the scattering
operator:

!equation id=discontinuous_sn_wf
\underbrace{\Big( \psi_{j},\, \frac{\partial}{\partial t}\frac{\Psi^{k}_{g,n}}{v_{g}} \Big)_{D^{k}}}_{\text{Kernel}}
+ \underbrace{\Big\langle \psi_{j},\, \hat{n}\cdot\hat{\Omega}_{n}\mathcal{F}^{*}_{g} \Big\rangle_{F^{k}}}_{\text{DG Kernel or Boundary Condition}}
- \underbrace{\Big( \vec{\nabla}\psi_{j}\cdot\hat{\Omega}_{n},\, \Psi^{k}_{g,n} \Big)_{D^{k}}}_{\text{Kernel}}\\
+ \underbrace{\Big( \psi_{j},\, \Sigma_{r,\,g}\Psi^{k}_{g,n} \Big)_{D^{k}}}_{\text{Kernel}}
- \underbrace{\Big( \psi_{j},\, \sum_{g' = 1}^{G}\Sigma_{s,\, g'\rightarrow g}\sum_{l = 0}^{L}\frac{2l + 1}{4\pi} f_{g'\rightarrow g,\, l}\sum_{m = -l}^{l}Y_{l,m}(\hat{\Omega}_{n})\Phi_{g',l,m} \Big)_{D^{k}}}_{\text{Kernel}}
- \underbrace{\Big( \psi_{j},\, S_{g,n} \Big)_{D^{k}}}_{\text{Kernel or Dirac Kernel}} = 0

The flux moments $\Phi_{g,l,m}$ are computed using an angular quadrature to treat the integral:

!equation id=angular_quad_final
\Phi_{g,l,m} = \sum_{i = 1}^{N} w_{n} Y_{l,m}(\hat{\Omega}_{n}) \Psi_{g,n}

### Neutron Transport Physics Objects

The the MooseObjects which implement each term in [!eqref](saaf_sn_wf),
[!eqref](discontinuous_sn_wf) and [!eqref](angular_quad_final) are summarized
below in [Table 2](#moose_objects).

!table id=moose_objects caption=MooseObjects corresponding to terms in [!eqref](saaf_sn_wf) and [!eqref](discontinuous_sn_wf).
| Object Name | Type | Implemented Term | Used If |
| - | - | - | - |
| [NeutronFluxMoment](source/auxkernels/NeutronFluxMoment.md) | AuxKernel | $\Phi_{g,l,m} = \sum_{i = 1}^{N} w_{n} Y_{l,m}(\hat{\Omega}_{n}) \Psi_{g,n}$ | SAAF or Upwinding |
| [ADSNMatSourceBC](source/bcs/ADSNMatSourceBC.md) | BC | $\Big\langle \psi_{j},\, \hat{n}\cdot\hat{\Omega}_{n}\Psi^{h}_{g,n} \Big\rangle_{\Gamma}$ or $\Big\langle \psi_{j},\, \hat{n}\cdot\hat{\Omega}_{n}\mathcal{F}^{*}_{g} \Big\rangle_{F^{k}}$ | SAAF or Upwinding |
| [ADSNReflectiveBC](source/bcs/ADSNReflectiveBC.md) | BC | $\Big\langle \psi_{j},\, \hat{n}\cdot\hat{\Omega}_{n}\Psi^{h}_{g,n} \Big\rangle_{\Gamma}$ or $\Big\langle \psi_{j},\, \hat{n}\cdot\hat{\Omega}_{n}\mathcal{F}^{*}_{g} \Big\rangle_{F^{k}}$ | SAAF or Upwinding |
| [ADSNVacuumBC](source/bcs/ADSNVacuumBC.md) | BC | $\Big\langle \psi_{j},\, \hat{n}\cdot\hat{\Omega}_{n}\Psi^{h}_{g,n} \Big\rangle_{\Gamma}$ or $\Big\langle \psi_{j},\, \hat{n}\cdot\hat{\Omega}_{n}\mathcal{F}^{*}_{g} \Big\rangle_{F^{k}}$ | SAAF or Upwinding |
| [ADSNRemoval](source/kernels/ADSNRemoval.md) | Kernel | $\Big( \psi_{j},\, \Sigma_{r,\,g}\Psi^{h}_{g,n} \Big)_{V}$ | SAAF or Upwinding |
| [ADSAAFTimeDerivative](source/kernels/ADSAAFTimeDerivative.md) | Kernel | $\Big( \psi_{j} + \tau_{g}\vec{\nabla}\psi_{j}\cdot\hat{\Omega}_{n},\, \frac{\partial}{\partial t}\frac{\Psi^{h}_{g,n}}{v_{g}} \Big)_{V}$ | SAAF\\ Transient |
| [ADSAAFStreaming](source/kernels/ADSAAFStreaming.md) | Kernel | $\Big( \vec{\nabla}\psi_{j}\cdot\hat{\Omega}_{n},\, \tau_{g}\hat{\Omega}_{n}\cdot\vec{\nabla}\Psi_{g,n}^{h} + (\tau_{g}\Sigma_{r,\, g} - 1)\Psi_{g,n}^{h} \Big)_{V}$ | SAAF |
| [ADSAAFScattering](source/kernels/ADSAAFScattering.md) | Kernel | $-\Big( \psi_{j} + \tau_{g}\vec{\nabla}\psi_{j}\cdot\hat{\Omega}_{n},\, \sum_{g' = 1}^{G}\Sigma_{s,\, g'\rightarrow g}\sum_{l = 0}^{L}\frac{2l + 1}{4\pi} f_{g'\rightarrow g,\, l}\sum_{m = -l}^{l}Y_{l,m}(\hat{\Omega}_{n})\Phi_{g',l,m} \Big)_{V}$ | SAAF\\ Assemble full scattering matrix |
| [ADSAAFExternalScattering](source/kernels/ADSAAFExternalScattering.md) | Kernel | $-\Big( \psi_{j} + \tau_{g}\vec{\nabla}\psi_{j}\cdot\hat{\Omega}_{n},\\ \sum_{g' = 1,g'\neq g}^{G}\Sigma_{s,\, g'\rightarrow g}\sum_{l = 0}^{L}\frac{2l + 1}{4\pi} f_{g'\rightarrow g,\, l}\sum_{m = -l}^{l}Y_{l,m}(\hat{\Omega}_{n})\Phi_{g',l,m} \Big)_{V}$ | SAAF\\ Source iteration |
| [ADSAAFInternalScattering](source/kernels/ADSAAFInternalScattering.md) | Kernel | $-\Big( \psi_{j} + \tau_{g}\vec{\nabla}\psi_{j}\cdot\hat{\Omega}_{n},\, \Sigma_{s,\, g\rightarrow g}\sum_{l = 0}^{L}\frac{2l + 1}{4\pi} f_{g\rightarrow g,\, l}\sum_{m = -l}^{l}Y_{l,m}(\hat{\Omega}_{n})\Phi_{g,l,m} \Big)_{V}$ | SAAF\\ Source iteration |
| [ADSAAFMaterialSource](source/kernels/ADSAAFMaterialSource.md) | Kernel | $- \Big( \psi_{j} + \tau_{g}\vec{\nabla}\psi_{j}\cdot\hat{\Omega}_{n},\, S_{g,n} \Big)_{V}$ | SAAF |
| [SAAFAnisoPointSource](source/dirackernels/SAAFAnisoPointSource.md) | Dirac Kernel | $- \Big( \psi_{j} + \tau_{g}\vec{\nabla}\psi_{j}\cdot\hat{\Omega}_{n},\, S_{g,n} \Big)_{V}$ | SAAF\\ Anisotropic Point Source |
| [SAAFIsoPointSource](source/dirackernels/SAAFIsoPointSource.md) | Dirac Kernel | $- \Big( \psi_{j} + \tau_{g}\vec{\nabla}\psi_{j}\cdot\hat{\Omega}_{n},\, S_{g,n} \Big)_{V}$ | SAAF\\ Isotropic Point Source |
| [ADDFEMTimeDerivative](source/kernels/ADDFEMTimeDerivative.md) | Kernel | $\Big( \psi_{j},\, \frac{\partial}{\partial t}\frac{\Psi^{k}_{g,n}}{v_{g}} \Big)_{D^{k}}$ | Upwinding |
| [ADDFEMStreaming](source/kernels/ADDFEMStreaming.md) | Kernel | $-\Big( \vec{\nabla}\psi_{j}\cdot\hat{\Omega}_{n},\, \Psi^{k}_{g,n} \Big)_{D^{k}}$ | Upwinding |
| [ADDFEMUpwinding](source/dgkernels/ADDFEMUpwinding.md) | DGKernel | $\Big\langle \psi_{j},\, \hat{n}\cdot\hat{\Omega}_{n}\mathcal{F}^{*}_{g} \Big\rangle_{F^{k}}$ | Upwinding |
| [ADDFEMScattering](source/kernels/ADDFEMScattering.md) | Kernel | $-\Big( \psi_{j},\, \sum_{g' = 1}^{G}\Sigma_{s,\, g'\rightarrow g}\sum_{l = 0}^{L}\frac{2l + 1}{4\pi} f_{g'\rightarrow g,\, l}\sum_{m = -l}^{l}Y_{l,m}(\hat{\Omega}_{n})\Phi_{g',l,m} \Big)_{D^{k}}$ | Upwinding\\ Assemble full scattering matrix |
| [ADDFEMExternalScattering](source/kernels/ADDFEMExternalScattering.md) | Kernel | $-\Big( \psi_{j},\, \sum_{g' = 1,g'\neq g}^{G}\Sigma_{s,\, g'\rightarrow g}\sum_{l = 0}^{L}\frac{2l + 1}{4\pi} f_{g'\rightarrow g,\, l}\sum_{m = -l}^{l}Y_{l,m}(\hat{\Omega}_{n})\Phi_{g',l,m} \Big)_{D^{k}}$ | Upwinding\\ Source iteration |
| [ADDFEMInternalScattering](source/kernels/ADDFEMInternalScattering.md) | Kernel | $-\Big( \psi_{j},\, \Sigma_{s,\, g\rightarrow g}\sum_{l = 0}^{L}\frac{2l + 1}{4\pi} f_{g\rightarrow g,\, l}\sum_{m = -l}^{l}Y_{l,m}(\hat{\Omega}_{n})\Phi_{g,l,m} \Big)_{D^{k}}$ | Upwinding\\ Source iteration |
| [ADDFEMMaterialSource](source/kernels/ADDFEMMaterialSource.md) | Kernel | $-\Big( \psi_{j},\, S_{g,n} \Big)_{D^{k}}$ | Upwinding |
| [DFEMAnisoPointSource](source/dirackernels/DFEMAnisoPointSource.md) | Dirac Kernel | $-\Big( \psi_{j},\, S_{g,n} \Big)_{D^{k}}$ | Upwinding\\ Anisotropic Point Source |
| [DFEMIsoPointSource](source/dirackernels/DFEMIsoPointSource.md) | Dirac Kernel | $-\Big( \psi_{j},\, S_{g,n} \Big)_{D^{k}}$ | Upwinding\\ Isotropic Point Source |

More information of the implementation of Moose Actions and Materials for Gnat
can be found in the
[transport system syntax page](syntax/NeutronActivationStudy/TransportSystem/index.md)
and the [material source pages](source/index.md).

### Neutron Transport Solver Objects

## Mass Transport
