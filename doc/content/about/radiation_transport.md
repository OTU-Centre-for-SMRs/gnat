## The Radiation Transport Equation

Gnat solves the multi-group Radiation Transport Equation (RTE) where the scattering kernel and external sources have been expanded in spherical harmonics:

!equation id=RTE
\frac{\partial}{\partial t}\frac{\Psi_{g}(\vec{r}, \hat{\Omega}, t)}{v_{g}} + \hat{\Omega}\cdot\vec{\nabla}\Psi_{g}(\vec{r}, \hat{\Omega}, t) + \Sigma_{t,\,g}(\vec{r})\Psi_{g}(\vec{r}, \hat{\Omega}, t) =
\sum_{g' = 1}^{G}\sum_{l = 0}^{L_{sc}}\frac{2l + 1}{4\pi}\Sigma_{s, g', l}(\vec{r}, t)\sum_{m = -l}^{l}Y_{l,m}(\hat{\Omega})\Phi_{g',l,m}(\vec{r}, t)
+ S(\vec{r},\hat{\Omega}, t).

In [!eqref](RTE) $\vec{r}$ is a position vector ($cm$), $\hat{\Omega}$ is the angular direction of travel ($st$), $t$ is time ($s$), $g$ is the radiation energy group, $G$ is the total number of energy groups, $\Psi_{g}(\vec{r}, \hat{\Omega}, t)$ is the angular radiation flux for group $g$ ($cm^{-2}s^{-1}st^{-1}$), $v_{g}$ is the radiation speed for a given energy group $g$ ($cm\, s^{-1}$), $\Sigma_{t,g}(\vec{r}, t)$ is the total macroscopic cross-section for group $g$ ($cm^{-1}$), $L_{sc}$ is the degree of the Legendre expansion on the scattering kernel, $\Sigma_{s, g', l}(\vec{r}, t)$ are the Legendre moments of the macroscopic scattering cross-section ($cm^{-1}st^{-1}$), $Y_{l,m}(\hat{\Omega})$ are the real spherical harmonics ($st^{-1}$), $\Phi_{g, l, m}(\vec{r}, t)$ are the spherical harmonic moments of the angular flux ($cm^{-2}s^{-1}$), $S_{g}(\vec{r}, \hat{\Omega}, t)$ is a group angular neutron source for group $g$ ($cm^{-3}s^{-1}st^{-1}$). $S_{g}(\vec{r}, \hat{\Omega}, t)$ has the following definition in Gnat:

!equation id=RTE_src
S_{g}(\vec{r},\hat{\Omega}, t)
= \frac{\chi_{g}(\vec{r})}{4\pi}\sum_{g'=1}^{G}\nu\Sigma_{f, g'}(\vec{r}, t)\Phi_{g', 0, 0}(\vec{r}, t)
+ \sum_{l = 0}^{L_{ext}}\frac{2l + 1}{4\pi}\sum_{m = -l}^{l}Y_{l,m}(\hat{\Omega})S_{g, l, m}^{ext}(\vec{r}, t),

where $\chi_{g}(\vec{r})$ is the fission spectra for group $g$, $\nu\Sigma_{f, g'}(\vec{r}, t)$ is the neutron production cross-section for group $g$ ($cm^{-1}$), $\Phi_{g', 0, 0}(\vec{r}, t)$ is the scalar neutron flux for group $g$ ($cm^{-2}s^{-1}$), and $L_{ext}$ is the degree of the spherical harmonics expansion of the group external source $S_{g, l, m}^{ext}(\vec{r}, t)$ ($cm^{-3}s^{-1}$). The fission source is only relevant to neutron transport problems, and is therefore excluded in gamma photon transport simulations.

The boundary conditions for [!eqref](RTE) are the following:

!equation id=RTE_bc
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

## Angular Discretization

## Spatial Discretization and Stabilization
