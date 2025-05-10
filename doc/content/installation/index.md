# Getting Started

The compilation instructions for Gnat are identical to that of any other MOOSE application,
with two different options for building and installing. The first option is to build using
the MOOSE conda environment, which contains pre-built binaries and libraries for all of MOOSE's
dependencies. This is the easiest, but the sparratic development of Gnat may result in
a mismatch between the most recent MOOSE conda environment and Gnat's MOOSE submodule
which often results in build errors. To avoid these errors we recommend building Gnat, MOOSE,
and MOOSE's dependencies from source instead of using the MOOSE conda environment.

# Build Instructions

1. [Using the MOOSE conda environment](moose_conda.md);
2. [Building entirely from source](from_source.md).
