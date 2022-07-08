# Governing Equations

## Neutron Transport

Gnat solves several governing equations. The first and chief among them is the
multi-group neutron transport equation (NTE) with the fission sources omitted.
The NTE has the following form:

!equation id=nte
\frac{\partial}{\partial t}\frac{\Psi_{g}(\vec{r}, \hat{\Omega}, t)}{v_{g}} + \hat{\Omega}\cdot\vec{\nabla}\Psi_{g}(\vec{r}, \hat{\Omega}, t) + \Sigma_{r,\,g}(\vec{r})\Psi_{g}(\vec{r}, \hat{\Omega}, t) =
\sum_{g' = 1}^{G}\int_{4\pi}\Sigma_{s,\, g'\rightarrow g}f_{g'\rightarrow g}(\vec{r}, \hat{\Omega}'\cdot\hat{\Omega})\Psi_{g'}(\vec{r}, \hat{\Omega}', t)\, d\hat{\Omega}'
+ S_{g}(\vec{r}, \hat{\Omega}, t)

where:

- $\vec{r}$ is a position vector ($cm$).
- $\hat{\Omega}$ is the angular direction of travel ($st$).
- $t$ is time ($s$).
- $g$ is the neutron energy group.
- $G$ is the total number of groups.
- $\Psi_{g}(\vec{r}, \hat{\Omega}, t)$ is the neutron angular flux for energy group $g$
  ($cm^{-2}s^{-1}st^{-1}$).
- $v_{g}$ is the neutron speed for a given energy group $g$ ($cm\, s^{-1}$).
- $\Sigma_{r,\,g}(\vec{r}, t)$ is the macroscopic neutron removal cross-section for group $g$ ($cm^{-1}$).
- $\Sigma_{s,\, g'\rightarrow g}f_{g'\rightarrow g}(\vec{r}, \Omega'\cdot\Omega, t)$
  is the scattering cross-section and phase function from group $g'$ into group $g$,
  and from angle $\hat{\Omega}'$ into angle $\hat{\Omega}$ \\ ($cm^{-1}st^{-1}$).
- $S_{g}(\vec{r}, \hat{\Omega}, t)$ is a group angular neutron source for group $g$
  ($cm^{-3}s^{-1}st^{-1}$).

The boundary conditions for [!eqref](nte) are the following:

!equation id=nte_bc
\Psi_{g}(\vec{r}, \hat{\Omega}, t) = \Psi_{inc,\, g}(\vec{r}, \hat{\Omega}, t) + \alpha_{s,\, g}(\vec{r})\Psi_{g}(\vec{r}, \hat{\Omega}_{r}, t),\\
\vec{r}\in\Gamma\text{ and } \hat{n}\cdot\hat{\Omega} < 0

where

- $\hat{\Omega}_{r}$ is the specular reflection direction ($st$). $\hat{\Omega}_{r} = \hat{\Omega} - (2\hat{\Omega}\cdot\hat{n})\hat{n}$
- $\Psi_{inc,\, g}(\vec{r}, \hat{\Omega}, t)$ is the incoming angular neutron flux for energy group $g$
  ($cm^{-2}s^{-1}st^{-1}$).
- $\alpha_{s,\, g}(\vec{r})$ is the specular reflection albedo (dimensionless).
- $\Psi_{g}(\vec{r}, \hat{\Omega}_{r}, t)$ is the specularly reflected angular neutron flux for energy group $g$
  ($cm^{-2}s^{-1}st^{-1}$).
- $\hat{n}$ is the boundary normal vector.
- $\Gamma$ is the boundary surface.

The discontinuous weak form of [!eqref](nte) is given by the following (with functional
notation omitted for brevity):

!equation id=nte_wf
\underbrace{\Big( \psi_{j},\, \frac{\partial}{\partial t}\frac{\Psi^{k}_{g}}{v_{g}} \Big)_{D^{k}}}_{\text{Time Kernel}}
+ \underbrace{\sum_{i} \Big\langle \psi_{j},\, \hat{n}\cdot\hat{\Omega}\Psi^{h}_{g}\Big|_{F^{k}_{\mathcal{I},\,i}-} \Big\rangle_{F^{k}_{\mathcal{I},\,i}}
+ \sum_{i} \Big\langle \psi_{j},\, \hat{n}\cdot\hat{\Omega}\Psi^{h}_{g}\Big|_{F^{k}_{\mathcal{O},\,i}+} \Big\rangle_{F^{k}_{\mathcal{O},\,i}}}_{\text{Discontinuous Galerkin Streaming Kernel}}
- \underbrace{\Big( \vec{\nabla}\psi_{j}\cdot\hat{\Omega},\, \Psi^{k}_{g} \Big)_{D^{k}}}_{\text{Streaming Kernel}}\\
+ \underbrace{\Big( \psi_{j},\, \Sigma_{r,\,g}\Psi^{k}_{g} \Big)_{D^{k}}}_{\text{Removal Kernel}}
- \underbrace{\Big( \psi_{j},\, \sum_{g' = 1}^{G}\int_{4\pi}\Sigma_{s,\, g'\rightarrow g}f_{g'\rightarrow g}\Psi^{k}_{g'}\, d\hat{\Omega}' \Big)_{D^{k}}}_{\text{Scattering Kernel}}
- \underbrace{\Big( \psi_{j},\, S_{g} \Big)_{D^{k}}}_{\text{Source Kernel}} = 0

The weak forms of [!eqref](nte_bc) are given by the following (with functional notation omitted for brevity):

!equation id=nte_bc_wf_r
\Big\langle \psi_{j},\, \hat{n}\cdot\hat{\Omega}\alpha_{s,\, g}\Psi_{r,\,g}^{k}\Big\rangle_{\Gamma_{r}} = 0,\quad\hat{n}\cdot\hat{\Omega} \leq 0
\\
\underbrace
{
\Big\langle \psi_{j},\, \hat{n}\cdot\hat{\Omega}\Psi_{g}^{k}\Big\rangle_{\Gamma_{r}} = 0,\quad\hat{n}\cdot\hat{\Omega} > 0
}_{\text{Reflective Boundary Condition}}

for the reflective boundary condition and

!equation id=nte_bc_wf_i
\Big\langle \psi_{j},\, \hat{n}\cdot\hat{\Omega}\Psi_{inc,\, g}\Big\rangle_{\Gamma_{i}} = 0,\quad\hat{n}\cdot\hat{\Omega} \leq 0\\
\underbrace{\Big\langle \psi_{j},\, \hat{n}\cdot\hat{\Omega}\Psi_{g}\Big\rangle_{\Gamma_{i}} = 0,\quad\hat{n}\cdot\hat{\Omega} > 0}_{\text{Incoming Flux Boundary Condition}}

for the incoming flux boundary condition. The combined boundary condition has been
decomposed into two boundary conditions for reflective boundaries $\Gamma_{r}$ and
incoming flux bondaries $\Gamma_{i}$ such that $\Gamma = \Gamma_{r}\cup\Gamma_{i}$.
A special case of an incoming flux boundary condition where $\Psi_{inc,\, g}(\vec{r}, \hat{\Omega}, t) = 0$
is the vacuum boundary condition.

The discrete ordiantes (+S@n@+) method is used to discretize the angular variable in
[!eqref](nte_wf), [!eqref](nte_bc_wf_r), and [!eqref](nte_bc_wf_i). The scattering kernel is expanded in Legendre
polynomials to yield a spherical harmonics representation for anisotropic scattering.
These methods are discussed in the [angular approach section](angular_approach.md).

## Mass Transport and Depletion
