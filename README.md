Gnat
=====

Gnat is an open-source radiation transport and mobile nuclide depletion solver built within the MOOSE
finite-element framework. Gnat has been designed for the simulation of neutron
activation of gases and aerosols, acting as the source-term solver for Caribou:
a multiphysics health physics and environmental impact code.

## Installation

Installation instructions can be found
[here](https://github.com/nuclearkevin/gnat/blob/master/doc/content/getting_started/installation.md).

## Documentation

Documentation of Gnat is still a work in progress. The current iteration of the docs
can be obtained by building the documentation website using the following series
of terminal commands (assuming the installation instructions have been followed):

```language=bash
mamba activate moose
cd gnat/doc
./moosedocs.py build --serve
```

Afterwards, navigate to [here](http://127.0.0.1:8000/source/index.html) with your browser of
choice to view the source file documentation. To view a breakdown of the governing
equations, visit [this page](http://127.0.0.1:8000/about/equations.html). The tutorials (when
written) can be found [here](http://127.0.0.1:8000/getting_started/tutorials.html).

The documentation website can be closed by entering `Ctrl + c` in the
terminal which built the documentation.

New users of Gnat and/or CARIBOU are encourages to visit the
[MOOSE website](https://mooseframework.inl.gov/) and the
[MOOSE tutorials](https://mooseframework.inl.gov/getting_started/examples_and_tutorials/index.html)
to learn about how MOOSE-based applications are structured and to get a feel for
using and developing a MOOSE-based application.

## Acknowledgements

We would like to acknowledge the support of the University Network of Excellence in Nuclear Engineering (UNENE) and the Natural Sciences and Engineering Research Council of Canada (NSERC) through the Industrial Research Chairs program [funding reference number IRCPJ 549979-19]. 

![NSERC Logo](https://www.nserc-crsng.gc.ca/img/logos/img-logo2-en.png)
