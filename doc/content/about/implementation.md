# Numerical Implementation in Moose

Moose is a Galerkin finite element framework which provides several plug-and-play
systems to assist in the implementation of partial differential equations. Gnat
takes advantage of the kernel, auxiliary kernel, Dirac kernel, boundary condition, material and action systems. The different
systems used by Gnat are summarized below:

!table id=used_moose_systems caption=Description of Moose systems used by Gnat.
| System Name | Description |
| - | - |
| Kernel | Evaluates the non-linear residual contribution of a term in the system of coupled partial differential equations for a single non-linear variable. |
| Auxiliary Kernel (AuxKernel) | Evaluates the value of a variable auxiliary to the problem.  |
| Dirac Kernel | Evaluates the non-linear residual contribution of a point source term in the system of coupled partial differential equations for a single non-linear variable. |
| Boundary Condition (BC) | Evaluates the non-linear residual contribution of the implicit boundary conditions in the system of coupled partial differential equations for a single non-linear variable. |
| Material | Provides domain-dependant material properties to other Moose systems. |
| Action | A system which allows Moose users to specialize and automate components of how Moose sets up, solves, and post-processes the problem.  |

## Neutral Particle Transport Implementation

Each residual contribution in a partial differential equation can be assigned to
a specific Moose system. The [SAAF weak form](stabilization.md#saaf) of the
Boltzmann transport equation with the
[+S@N@+ disretization](nte_angular_approach.md#sn_disc) applied in angle and
[+P@N@+ treatment](nte_angular_approach.md#scattering) applied to the scattering
operator is:

!equation id=saaf_sn_wf
\underbrace{\Big( \psi_{j} + \tau_{g}\vec{\nabla}\phi_{j}\cdot\hat{\Omega}_{n},\, \frac{\partial}{\partial t}\frac{\Psi^{h}_{g,n}}{v_{g}} \Big)_{V}}_{\text{Kernel}}
+ \underbrace{\Big\langle \psi_{j},\, \hat{n}\cdot\hat{\Omega}_{n}\Psi^{h}_{g,n} \Big\rangle_{\Gamma}}_{\text{Boundary Condition}}
+ \underbrace{\Big( \vec{\nabla}\psi_{j}\cdot\hat{\Omega}_{n},\, \tau_{g}\hat{\Omega}_{n}\cdot\vec{\nabla}\Psi_{g,n}^{h} + (\tau_{g}\Sigma_{r,\, g} - 1)\Psi_{g,n}^{h} \Big)_{V}}_{\text{Kernel}}\\
+ \underbrace{\Big( \psi_{j},\, \Sigma_{r,\,g}\Psi^{h}_{g,n} \Big)_{V}}_{\text{Kernel}}
- \underbrace{\Big( \psi_{j} + \tau_{g}\vec{\nabla}\phi_{j}\cdot\hat{\Omega}_{n},\, \sum_{g' = 1}^{G}\Sigma_{s,\, g'\rightarrow g}\sum_{l = 0}^{L}\frac{2l + 1}{4\pi} f_{g'\rightarrow g,\, l}\sum_{m = -l}^{l}Y_{l,m}(\hat{\Omega}_{n})\Phi_{g',l,m} \Big)_{V}}_{\text{Kernel}}\\
-\underbrace{\Big( \psi_{j} + \tau_{g}\vec{\nabla}\psi_{j}\cdot\hat{\Omega}_{n}, \frac{\chi_{g}}{4\pi}\sum_{g' = 1}^{G}\nu\Sigma_{f,\, g'}\Phi_{g',0,0} \Big)_{V}}_{\text{Kernel}}
- \underbrace{\Big( \psi_{j} + \tau_{g}\vec{\nabla}\phi_{j}\cdot\hat{\Omega}_{n},\, S_{g,n} \Big)_{V}}_{\text{Kernel or Dirac Kernel}} = 0

The flux moments $\Phi_{g,l,m}$ are computed using an angular quadrature to treat the integral:

!equation id=angular_quad_final
\Phi_{g,l,m} = \sum_{i = 1}^{N} w_{n} Y_{l,m}(\hat{\Omega}_{n}) \Psi_{g,n}

### Neutral Particle Transport Physics Objects

The the MooseObjects which implement each term in [!eqref](saaf_sn_wf) and [!eqref](angular_quad_final) are summarized
below in [Table 2](#moose_objects).

!table id=moose_objects caption=MooseObjects corresponding to terms in [!eqref](saaf_sn_wf).
| Object Name | Type | Implemented Term | Used If |
| - | - | - | - |
| [ParticleFluxMoment](source/auxkernels/ParticleFluxMoment.md) | AuxKernel | $\Phi_{g,l,m} = \sum_{i = 1}^{N} w_{n} Y_{l,m}(\hat{\Omega}_{n}) \Psi_{g,n}$ | SAAF |
| [SNSourceBC](source/bcs/SNSourceBC.md) | BC | $\Big\langle \psi_{j},\, \hat{n}\cdot\hat{\Omega}_{n}\Psi^{h}_{g,n} \Big\rangle_{\Gamma}$ | SAAF\\ Fixed source boundary |
| [SNReflectiveBC](source/bcs/SNReflectiveBC.md) | BC | $\Big\langle \psi_{j},\, \hat{n}\cdot\hat{\Omega}_{n}\Psi^{h}_{g,n} \Big\rangle_{\Gamma}$ | SAAF\\ Reflective boundary |
| [SNVacuumBC](source/bcs/SNVacuumBC.md) | BC | $\Big\langle \psi_{j},\, \hat{n}\cdot\hat{\Omega}_{n}\Psi^{h}_{g,n} \Big\rangle_{\Gamma}$ | SAAF\\ Non-reentrant boundary |
| [SNRemoval](source/kernels/SNRemoval.md) | Kernel | $\Big( \psi_{j},\, \Sigma_{r,\,g}\Psi^{h}_{g,n} \Big)_{V}$ | SAAF |
| [SAAFTimeDerivative](source/kernels/SAAFTimeDerivative.md) | Kernel | $\Big( \psi_{j} + \tau_{g}\vec{\nabla}\psi_{j}\cdot\hat{\Omega}_{n},\, \frac{\partial}{\partial t}\frac{\Psi^{h}_{g,n}}{v_{g}} \Big)_{V}$ | SAAF\\ Transient |
| [SAAFStreaming](source/kernels/SAAFStreaming.md) | Kernel | $\Big( \vec{\nabla}\psi_{j}\cdot\hat{\Omega}_{n},\, \tau_{g}\hat{\Omega}_{n}\cdot\vec{\nabla}\Psi_{g,n}^{h} + (\tau_{g}\Sigma_{r,\, g} - 1)\Psi_{g,n}^{h} \Big)_{V}$ | SAAF |
| [SAAFScattering](source/kernels/SAAFScattering.md) | Kernel | $-\Big( \psi_{j} + \tau_{g}\vec{\nabla}\psi_{j}\cdot\hat{\Omega}_{n},\,\\ \sum_{g' = 1}^{G}\Sigma_{s,\, g'\rightarrow g}\sum_{l = 0}^{L}\frac{2l + 1}{4\pi} f_{g'\rightarrow g,\, l}\sum_{m = -l}^{l}Y_{l,m}(\hat{\Omega}_{n})\Phi_{g',l,m} \Big)_{V}$ | SAAF\\ User enables off-diagonal jacobian contributions |
| [SAAFMomentScattering](source/kernels/SAAFMomentScattering.md) | Kernel | $-\Big( \psi_{j} + \tau_{g}\vec{\nabla}\psi_{j}\cdot\hat{\Omega}_{n},\\ \sum_{g' = 1}^{G}\Sigma_{s,\, g'\rightarrow g}\sum_{l = 0}^{L}\frac{2l + 1}{4\pi} f_{g'\rightarrow g,\, l}\sum_{m = -l}^{l}Y_{l,m}(\hat{\Omega}_{n})\Phi_{g',l,m} \Big)_{V}$ | SAAF\\ User disables off-diagonal jacobian contributions |
| [SAAFMomentFission](source/kernels/SAAFMomentFission.md) | Kernel | $-\Big( \psi_{j} + \tau_{g}\vec{\nabla}\psi_{j}\cdot\hat{\Omega}_{n}, \frac{\chi_{g}}{4\pi}\sum_{g' = 1}^{G}\nu\Sigma_{f,\, g'}\Phi_{g',0,0} \Big)_{V}$ | SAAF\\ Particle is neutron\\ User enables fission\\ User disables off-diagonal jacobian contributions |
| [SAAFVolumeSource](source/kernels/SAAFVolumeSource.md) | Kernel | $- \Big( \psi_{j} + \tau_{g}\vec{\nabla}\psi_{j}\cdot\hat{\Omega}_{n},\, S_{g,n} \Big)_{V}$ | SAAF\\ Block contains a constant volumetric source |
| [SAAFFieldSource](source/kernels/SAAFFieldSource.md) | Kernel | $- \Big( \psi_{j} + \tau_{g}\vec{\nabla}\psi_{j}\cdot\hat{\Omega}_{n},\, S_{g,n}(\vec{r}, t) \Big)_{V}$ | SAAF\\ Block contains a particle source that is a function of another variable |
| [SAAFPointSource](source/dirackernels/SAAFPointSource.md) | Dirac Kernel | $- \Big( \psi_{j} + \tau_{g}\vec{\nabla}\psi_{j}\cdot\hat{\Omega}_{n},\, \delta(\vec{r} - \vec{r}_{s})S_{g,n} \Big)_{V}$ | SAAF\\ Problem contains a constant point source |

More information of the implementation of Moose Actions and Materials for Gnat
can be found in the
[transport system syntax page](syntax/TransportSystems/index.md)
and the [material source pages](source/index.md).

### Neutral Particle Transport Solver Objects

## Trace Species Transport
