# Angular Discretization

## Discrete Ordinates

The weak form of the neutron transport equation ([!eqref](equations.md#nte)) is
discretized in the angular dimension with the use of the discrete ordinates (+S@n@+)
method. The NTE is solved for a collection of discrete ordinates (directions) $\hat{\Omega}_{n}$,
and a suitable integration method is used to compute the angular integral using
the solutions along the ordinates. The +S@n@+ weak form of the NTE
is given with the following:

!equation id=nte_wf_sn
\underbrace{\Big( \psi_{j},\, \frac{\partial}{\partial t}\frac{\Psi^{k}_{g,n}}{v_{g}} \Big)_{D^{k}}}_{\text{Time Kernel}}
+ \underbrace{+ \Big\langle \phi_{j},\, \hat{n}\cdot\hat{\Omega}\Psi^{h}_{g,n} \Big\rangle_{\Gamma}}_{\text{Implicit Boundary Conditions}}
- \underbrace{\Big( \vec{\nabla}\psi_{j}\cdot\hat{\Omega},\, \Psi^{k}_{g,n} \Big)_{D^{k}}}_{\text{Streaming Kernel}}\\
+ \underbrace{\Big( \psi_{j},\, \Sigma_{r,\,g}\Psi^{k}_{g,n} \Big)_{D^{k}}}_{\text{Removal Kernel}}
- \underbrace{\Big( \psi_{j},\, \sum_{g' = 1}^{G}\int_{4\pi}\Sigma_{s,\, g'\rightarrow g}f_{g'\rightarrow g}\Psi^{k}_{g'}\, d\hat{\Omega}' \Big)_{D^{k}}}_{\text{Scattering Kernel}}
- \underbrace{\Big( \psi_{j},\, S_{g, n} \Big)_{D^{k}}}_{\text{Source Kernel}} = 0

where:

- $\Psi^{k}_{g,n} = \Psi^{k}_{g}(\vec{r}, \hat{\Omega}_{n}, t)$.
- $S_{g, n}(\vec{r}, t) = S_{g}(\vec{r}, \hat{\Omega}_{n}, t)$.

## Scattering Treatment

The scattering cross-section and phase function can be expanded in Legendre polynomials $P_{l}$:

!equation
f_{g'\rightarrow g}(\vec{r}, \hat{\Omega}'\cdot\hat{\Omega}) \approx  
\sum_{l = 0}^{L}(2l + 1)f_{g'\rightarrow g,\, l}(\vec{r})
P_{l}(\hat{\Omega}'\cdot\hat{\Omega})

!equation
f_{g'\rightarrow g,\, l}(\vec{r}) = \int_{-1}^{1}f_{g'\rightarrow g}(\vec{r}, \mu_{0}) P_{l}(\mu_{0})\, d\mu_{0}

where $\mu_{0} = \hat{\Omega}'\cdot\hat{\Omega}$. This results in the modified scattering kernel as shown below:

!equation
\sum_{g' = 1}^{G}\int_{4\pi}\Sigma_{s,\, g'\rightarrow g}f_{g'\rightarrow g}\Psi_{g'}\, d\hat{\Omega}' \approx
\sum_{g' = 1}^{G}\Sigma_{s,\, g'\rightarrow g}(\vec{r})\sum_{l = 0}^{L}(2l + 1) f_{g'\rightarrow g,\, l}(\vec{r})
\int_{4\pi}P_{l}(\hat{\Omega}'\cdot\hat{\Omega})\Psi_{g'}(\vec{r}, \hat{\Omega}', t)\, d\hat{\Omega}'

The spherical harmonics addition theorem can be applied to remove the Legendre polynomial
$P_{l}(\hat{\Omega}'\cdot\hat{\Omega})$:

!equation id=scat_exp
\sum_{g' = 1}^{G}\int_{4\pi}\Sigma_{s,\, g'\rightarrow g}f_{g'\rightarrow g}\Psi_{g'}\, d\hat{\Omega}' \approx
\sum_{g' = 1}^{G}\Sigma_{s,\, g'\rightarrow g}(\vec{r})\sum_{l = 0}^{L}\frac{2l + 1}{4\pi} f_{g'\rightarrow g,\, l}(\vec{r})\sum_{m = -l}^{l}Y_{l,m}(\hat{\Omega})\Phi_{g',l,m}(\vec{r}, t)

!equation
\Phi_{g,l,m}(\vec{r}, t) = \int_{4\pi} Y_{l,m}(\hat{\Omega}') \Psi_{g}(\vec{r}, \hat{\Omega}', t)\, d\hat{\Omega}'

where:

- $l$ and $L$ are the degree and maximum degree of the Legendre expansion.
- $m$ is the order of the spherical harmonics expansion.
- $Y_{l,m}(\hat{\Omega})$ are the real spherical harmonics.
- $f_{g'\rightarrow g,\, l}(\vec{r})$ are the moments of the group-to-group
  scattering phase function.
- $\Phi_{g,l,m}(\vec{r}, t)$ are the moments of the group angular flux.

Substituting [!eqref](scat_exp) into the scattering kernel yields the following:

!equation id=ani_scat_kernel
\Big( \psi_{j},\, \sum_{g' = 1}^{G}\int_{4\pi}\Sigma_{s,\, g'\rightarrow g}f_{g'\rightarrow g}\Psi^{k}_{g'}\, d\hat{\Omega}' \Big)_{D^{k}} \approx
\underbrace{\Big( \psi_{j},\, \sum_{g' = 1}^{G}\Sigma_{s,\, g'\rightarrow g}\sum_{l = 0}^{L}\frac{2l + 1}{4\pi} f_{g'\rightarrow g,\, l}\sum_{m = -l}^{l}Y_{l,m}(\hat{\Omega}_{n})\Phi_{g',l,m} \Big)_{D^{k}}}_{\text{Scattering Kernel}}

In general the degree of the
Legendre expansion corresponds to the degree of anisotropy in the medium.

## Evaluating Flux Moments

The scattering expansion discussed above hinges on the evaluation of flux moments:

!equation
\Phi_{g,l,m}(\vec{r}, t) = \int_{4\pi} Y_{l,m}(\hat{\Omega}') \Psi_{g}(\vec{r}, \hat{\Omega}', t)\, d\hat{\Omega}'

The integral over $4\pi\,st$ of the unit sphere (all directions) is equivalent to the spherical integral:

!equation
\int_{4\pi} f(\hat{\Omega})\,d\hat{\Omega} = \int_{-1}^{1}\int_{0}^{2\pi}f(\mu,\omega)\,d\omega d\mu

For functions symmetrical about $\omega = \pi$ it can be shown that this integral is equal to:

!equation
\int_{4\pi} f(\hat{\Omega})\,d\hat{\Omega} = 2\int_{-1}^{1}\int_{-1}^{1}\frac{f(\hat{\Omega}(\mu,y))}{\sqrt{1-y^{2}}}\,dy d\mu

where the angular direction $\Omega$ is a function of $\mu$ and $y$. For functions
which are not symmetrical about $\omega = \pi$ we can decompose the integral
about $\pi$ such that both halves of the unit spheres are integrated completely:

!equation
\int_{4\pi} f(\hat{\Omega})\,d\hat{\Omega} = \int_{-1}^{1}\int_{-1}^{1}\frac{f(\hat{\Omega}(\mu,y))}{\sqrt{1-y^{2}}}\,dy d\mu + \int_{-1}^{1}\int_{-1}^{1}\frac{f(-\hat{\Omega}(\mu,y))}{\sqrt{1-y^{2}}}\,dy d\mu

This form allows us to evaluate angular integrals with a Gauss-Chebyshev product
quadrature, where the polar integral is approximated with a Gauss-Legendre quadrature
and the azimuthal integral is approximated with a Gauss-Chebyshev quadrature. This takes the form of:

!equation id=quad
\int_{-1}^{1}\int_{-1}^{1}\frac{f(\hat{\Omega}(\mu,y))}{\sqrt{1-y^{2}}}\,dy d\mu \approx
\sum_{i = 1}^{n_{L}}\sum_{j = 1}^{n_{C}} w_{i,\,L}w_{j,\,C} f(\hat{\Omega}(\mu_{i,\,L},y_{j,\,C}))

where:

- $n_{L}$ is the number of Gauss-Legendre quadrature points.
- $w_{i,\, L}$ is the $i$'th Gauss-Legendre weight.
- $\mu_{i,\, L}$ is the $i$'th Gauss-Legendre point.
- $n_{C}$ is the number of Gauss-Chebyshev quadrature points.
- $w_{j,\, C}$ is the $j$'th Gauss-Chebyshev weight.
- $y_{j,\, C}$ is the $j$'th Gauss-Chebyshev point.

$\mu_{i,\, L}$ are the roots of the $n_{L}$ degree Legendre polynomial, and
the weights $w_{i,\,L}$ are defined as:

!equation id=gl_quad
w_{i,\,L} = \frac{2}{(1-\mu_{i,\,L}^{2})(P'_{n_{L}}(\mu_{i,\,L}))}

$y_{j,\, C}$ and $w_{j,\, C}$ are:

!equation id=gc_quad
w_{j,\, C} = \frac{\pi}{n_{C}}\\
\,\\
y_{j,\, C} = \cos{(\omega_{j,\,C})},\, \omega_{j,\,C} = \frac{(2j - 1)\pi}{2n_{C}},\, j\in \left[ 1,\, n_{C} \right]

The order of the quadrature is $2n_{L}n_{C}$, as the angular flux is not guaranteed
to be symmetrical about $\omega = \pi$ and must therefore be evaluated twice for
positive and negative angular directions. We can chose the discrete ordinates
$\hat{\Omega}_{n}$ such that:

!equation
\hat{\Omega}_{n} = \hat{\Omega}(\mu_{i,\, L}, y_{j,\, C})

Allowing for the use of the discrete ordinates solutions to evaluate the required
flux moments.
