!config navigation breadcrumbs=False scrollspy=False

# HOME style=visibility:hidden;

!media gnat_logo.png
  style=width:75%;margin-left:auto;margin-right:auto;halign:center
  id=gnat_logo
  prefix=

# An open-source neutral particle transport and fluid depletion solver built on MOOSE
  class=center
  style=font-size:200%

Gnat provides an S@N@ radiation transport solver which is geared towards fluid depletion and field-based
radiation sources. This capability has been primarily designed with source-driven problems in mind
such as ex-core neutron activation, radiation shielding problems, and gamma photon transport from plumes.
However, neutron-induced fission reactions and eigenvalue solutions to the neutron transport equation (k@eff@)
are provided from the sake of completeness. The fluid activation and depletion models provide both finite element and finite
volume discretization schemes geared towards medium-fidelity simulations of dispersion-depletion problems.
As Gnat is built on MOOSE, it is capable of coupling with other MOOSE-based applications (such as
[Cardinal](https://cardinal.cels.anl.gov/index.html)) for higher-resolution simulations.

!col! small=12 medium=4 large=4 icon=get_app

### Getting Started class=center style=font-weight:200

!style halign=center
Gnat is open-source and the repository is available on [GitHub](https://github.com/nuclearkevin/gnat) -
learn how to [get started](getting_started/installation.md).
!col-end!

!col! small=12 medium=4 large=4 icon=group

### Tutorials and Documentation class=center style=font-weight:200

!style halign=center
Explore the [tutorials](tutorials/index.md) and [source documentation](source/index.md)
to learn about the capabilities of Gnat.
!col-end!

!col! small=12 medium=4 large=4 icon=settings

### Theory and Implementation class=center style=font-weight:200

!style halign=center
Review the [radiation transport](about/radiation_transport.md) and
[fluid depletion](about/mobile_depletion.md) models included in Gnat.

!col-end!

# Sample Uses
  class=center
  style=font-size:200%


!gallery! large=4
!card media/landing_page/subcritical_g25.png title=Criticality Calculations
Simulation of a graphite moderated natural uranium fueled subcritical assembly using the S@N@ radiation transport solver running in eigenvalue mode. Figure shows the thermal neutron flux.

!card media/landing_page/bwr_g2.png title=Ex-Core Radiation Transport
Simulation of the neutron fields in a concrete-walled containment system caused by a small BWR using the S@N@ radiation transport solver driven by a surface source. Figure shows the fast (group 2) scalar neutron flux.

!card media/landing_page/bwr_ar41.png title=Neutron Activation of Fluids
Transient activation of Ar-40 to Ar-41 in air near a small BWR using the neutron transport and fluid activation solvers. Figure shows the activity of Ar-41 in the containment system.
!gallery-end!

!gallery! large=6
!card media/landing_page/plume_cs137.png title=Radionuclide Plumes
Simulation of short length and time scale plume release from nuclear facilities with the tracer transport model. Figure shows the concentration of Cs-137 after 7.5 minutes of the release.

!card media/landing_page/plume_flux_log.png title=External Photon Dosimetry
Simulation of skyshine from the decay of Ba-137m in a plume release using the coupled tracer transport and photon transport solvers. Figure shows the monochromatic 662 keV photon flux.
!gallery-end!

# Acknowledgements
  class=center
  style=font-size:200%

The development of Gnat was supported by an Industrial Research Chair (IRC) program. We would
like to acknowledge the University Network of Excellence in Nuclear Engineering
(UNENE) and the Natural Sciences and Engineering Research Council of Canada (NSERC) for funding
this IRC position \[funding reference number IRCPJ 549979-19\].

!row! style=display:inline-flex;
!col! small=12 medium=4 large=3

!media otu_logo.png style=height:100%;display:block;halign:center;

!col-end!

!col! small=12 medium=4 large=3

!media nserc_symbol.png style=height:100%;padding:0px 0px;display:block;halign:center;

!col-end!

!col! small=12 medium=4 large=3

!media UNENE.png style=width:100%;display:block;halign:center;

!col-end!
!row-end!
