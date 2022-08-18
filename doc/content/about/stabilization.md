# Stabilization for the Finite Element Method

Both the neutron transport equation and the mass transport equations are
hyperbolic and are numerically unstable when solved using the finite element
method. There are several approaches to mitigate the oscillations that occur in
the neutron transport equation and mass transport problems dominated by
streaming and convection (respectively). The following sections discuss the approaches taken by
Gnat for both sets of governing equations, starting with the neutron transport
equation and ending with the convection-diffusion equations.

## Neutron Transport

Gnat provides two stabilization schemes for the neutron transport equation, both being described below.

### 1. The Self-Adjoint Angular Flux Method

The first approached implemented is the self-adjoint angular flux (SAAF) method
derived by (TODO: REFERENCE HERE). We can write the transport equation where
the collision term is on the left and all remaining terms are on the right:

!equation
\Sigma_{r,\,g}\Psi_{g} =
\sum_{g' = 1}^{G}\int_{4\pi}\Sigma_{s,\, g'\rightarrow g}f_{g'\rightarrow g}\Psi_{g'}\, d\hat{\Omega}'
+ S_{g}
- \frac{\partial}{\partial t}\frac{\Psi_{g}}{v_{g}}
- \hat{\Omega}\cdot\vec{\nabla}\Psi_{g}

Inverting the collision operator yields the following equation, which is known
as the angular flux equation:

!equation id=afe
\Psi_{g} = \frac{1}{\Sigma_{r,\,g}}
\Big(\sum_{g' = 1}^{G}\int_{4\pi}\Sigma_{s,\, g'\rightarrow g}f_{g'\rightarrow g}\Psi_{g'}\, d\hat{\Omega}'
+ S_{g}
- \frac{\partial}{\partial t}\frac{\Psi_{g}}{v_{g}}
- \hat{\Omega}\cdot\vec{\nabla}\Psi_{g}\Big)

[!eqref](afe) can be substituted into [!eqref](equations.md#nte_wf) and
algebraically simplified to yield the SAAF weak form of the neutron transport equation:

!equation id=saaf_wf_no_void
\Big( \phi_{j} + \frac{1}{\Sigma_{r,\,g}}\vec{\nabla}\phi_{j}\cdot\hat{\Omega},\, \frac{\partial}{\partial t}\frac{\Psi^{h}_{g}}{v_{g}} \Big)_{V}
+ \Big\langle \phi_{j},\, \hat{n}\cdot\hat{\Omega}\Psi^{h}_{g} \Big\rangle_{\Gamma}
+ \Big( \vec{\nabla}\phi_{j}\cdot\hat{\Omega},\, \frac{1}{\Sigma_{r,\,g}}\hat{\Omega}\cdot\vec{\nabla}\Psi_{g}^{h}\Big)_{V}
+ \Big( \phi_{j},\, \Sigma_{r,\,g}\Psi^{h}_{g} \Big)_{V}\\
- \Big( \phi_{j} + \frac{1}{\Sigma_{r,\,g}}\vec{\nabla}\phi_{j}\cdot\hat{\Omega},\, \sum_{g' 1}^{G}\int_{\mathcal{S}}\Sigma_{s,\, g'\rightarrow g}f_{g'\rightarrow g}\Psi^{h}_{g'}\, d\Omega' \Big)_{V}
- \Big( \phi_{j} + \frac{1}{\Sigma_{r,\,g}}\vec{\nabla}\phi_{j}\cdot\hat{\Omega},\, S_{g} \Big)_{V} = 0

[!eqref](saaf_wf_no_void) is second order and therefore stable for
continuous finite elements. However, the division by the removal cross-section
makes [!eqref](saaf_wf_no_void) ill-suited for regions where the removal
cross-section is void or close to being a void. To remedy this an alternative form of the
angular flux equation is proposed by (TODO: REFERENCE HERE):

!equation id=afe_void
\Psi_{g} = \Big(1 - \tau_{g}\Sigma_{r,\, g}\Big)\Psi_{g} +
\tau_{g}\Big(\sum_{g' = 1}^{G}\int_{\mathcal{S}}\Sigma_{s,\, g'\rightarrow g}f_{g'\rightarrow g}\Psi_{g'}\, d\Omega'
+ S_{g}
- \frac{\partial}{\partial t}\frac{\Psi_{g}}{v_{g}}
- \hat{\Omega}\cdot\vec{\nabla}\Psi_{g}\Big)

where $\tau_{g}$ is defined as:

!equation id=tau_g
\tau_{g} =
\begin{cases}
\frac{1}{c\Sigma_{r,\,g}}, & ch\Sigma_{r,\,g} \geq \zeta\\
\frac{h}{\zeta}, & ch\Sigma_{r,\,g} < \zeta
\end{cases}

h is the maximum vertex separation within a single finite element, c and $\zeta$
are constants to determine the maximum stabilization in a void region. Selecting
$\zeta = 0$ collapses [!eqref](afe_void) into [!eqref](afe). Substituting
[!eqref](afe_void) into [!eqref](equations.md#nte_wf) yields the void-treated SAAF weak form which is implemented in Gnat as the SAAF scheme:

!equation id=saaf_wf_void
\Big( \phi_{j} + \tau_{g}\vec{\nabla}\phi_{j}\cdot\hat{\Omega},\, \frac{\partial}{\partial t}\frac{\Psi^{h}_{g}}{v_{g}} \Big)_{V}
+ \Big\langle \phi_{j},\, \hat{n}\cdot\hat{\Omega}\Psi^{h}_{g} \Big\rangle_{\Gamma}
+ \Big( \vec{\nabla}\phi_{j}\cdot\hat{\Omega},\, \tau_{g}\hat{\Omega}\cdot\vec{\nabla}\Psi_{g}^{h} + (\tau_{g}\Sigma_{r,\, g} - 1)\Psi_{g}^{h} \Big)_{V}\\
+ \Big( \phi_{j},\, \Sigma_{r,\,g}\Psi^{h}_{g} \Big)_{V}
- \Big( \phi_{j} + \tau_{g}\vec{\nabla}\phi_{j}\cdot\hat{\Omega},\, \sum_{g' = 1}^{G}\int_{\mathcal{S}}\Sigma_{s,\, g'\rightarrow g}f_{g'\rightarrow g}\Psi^{h}_{g'}\, d\Omega' \Big)_{V}
- \Big( \phi_{j} + \tau_{g}\vec{\nabla}\phi_{j}\cdot\hat{\Omega},\, S_{g} \Big)_{V} = 0

### 2. The Upwinding Method

The upwinding approach for the neutron transport equation follows from the
discontinuous weak form of the neutron transport equation:

!equation id=discontinuous_wf
\Big( \phi_{j},\, \frac{\partial}{\partial t}\frac{\Psi^{k}_{g}}{v_{g}} \Big)_{D^{k}}
+ \Big\langle \phi_{j},\, \hat{n}\cdot\hat{\Omega}\mathcal{F}^{*}_{g} \Big\rangle_{F^{k}}
- \Big( \vec{\nabla}\phi_{j}\cdot\hat{\Omega},\, \Psi^{k}_{g} \Big)_{D^{k}}\\
+ \Big( \phi_{j},\, \Sigma_{r,\,g}\Psi^{k}_{g} \Big)_{D^{k}}
- \Big( \phi_{j},\, \sum_{g' = 1}^{G}\int_{\mathcal{S}}\Sigma_{s,\, g'\rightarrow g}f_{g'\rightarrow g}\Psi^{k}_{g'}\,d\Omega' \Big)_{D^{k}}
- \Big( \phi_{j},\, S_{g} \Big)_{D^{k}} = 0

where:

- $\Psi^{k}_{g}$ are the trial functions defined on the single finite element $k$.
  It is important to note that the trial functions are discontinuous over the full domain.
- $D^{k}$ is the volume of the finite element $k$ such that $D^{k}\in V$.
- $F^{k}$ is the bounding surface of the finite element $k$.
- $\mathcal{F}^{*}_{g}(\vec{r}) = \Psi^{k}_{g}(\vec{r}),\, \vec{r}\in F^{k}$
  is the interface numerical angular flux.

The upwinding scheme is comparatively simple to the SAAF scheme. For each finite
element $k$ the faces which compose the bounding surface $F^{k}$ are traversed.
Each face is classified as either upwind or downwind according to the
following scheme:

- Upwind faces are faces where $\hat{n}\cdot\hat{\Omega} \geq 0$.
- Downwpwind faces are faces where $\hat{n}\cdot\hat{\Omega} \lt 0$.

The numerical flux is then equal to the value of the neighbour's trial function
at the interface for downwind faces, and the value of the current element's
trial function at the interface.

## Mass Transport
