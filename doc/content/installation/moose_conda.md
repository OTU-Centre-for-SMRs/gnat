# Building with the MOOSE Conda Environment

## 1. Installing the MOOSE Environment

Gnat is a MOOSE-based application, and therefore requires an appropriate MOOSE
environment. More detailed instructions for installing Mambaforge3 (a dependency
manager) can be found on the MOOSE [installation page](https://mooseframework.inl.gov/getting_started/installation/conda.html). You can install Mambaforge3 with the following
shell commands:

Linux:

```language=bash
curl -L -O https://github.com/conda-forge/miniforge/releases/latest/download/Mambaforge-Linux-x86_64.sh
bash Mambaforge-Linux-x86_64.sh -b -p ~/mambaforge3
```

Macintosh without Apple Silicon:

```language=bash
curl -L -O https://github.com/conda-forge/miniforge/releases/latest/download/Mambaforge-MacOSX-x86_64.sh
bash Mambaforge-MacOSX-x86_64.sh -b -p ~/mambaforge3
```

Macintosh with Apple Silicon:

```language=bash
curl -L -O https://github.com/conda-forge/miniforge/releases/latest/download/Mambaforge-MacOSX-arm64.sh
bash Mambaforge-MacOSX-arm64.sh -b -p ~/mambaforge3
```

Once Mambaforge3 has been installed for your platform of choice, you need to
configure it to add the MOOSE environment. To do so, add Mambaforge3 to your system PATH:

```language=bash
export PATH=$HOME/mambaforge3/bin:$PATH
```

Then, add the INL public channel:

```language=bash
conda config --add channels https://conda.software.inl.gov/public
```

Finally, initialize the Mambaforge3 shell:

```language=bash
mamba init
```

and create+install the MOOSE environment:

```language=bash
conda create -n moose moose-dev=2025.05.13=mpich
conda activate moose
```

This will take some time, so feel free to walk away and get yourself a coffee (or
other beverage of your choice). Once it's complete deactivate your MOOSE environment
to ensure the installation is applied.

```language=bash
mamba deactivate
```

After the installation is finished you can move on to cloning Gnat and MOOSE.

## 2. Cloning Gnat and MOOSE

Gnat comes with the MOOSE framework included as a Git submodule. For users who wish
to develop applications alongside Gnat, use other sub-applications of Caribou, or
use other MOOSE applications: it's advised that you install MOOSE separately. This
can be done by following all of the detailed instructions on the MOOSE
[installation page](https://mooseframework.inl.gov/getting_started/installation/conda.html).

For users who only wish to use Gnat:

```language=bash
git clone https://github.com/OTU-Center-for-SMRs/gnat.git
cd gnat
git submodule update --init moose
```

For users who want to develop/use other MOOSE applications:

```language=bash
git clone https://github.com/OTU-Center-for-SMRs/gnat.git
```

Before continuing with the compilation process, ensure that the MOOSE environment
has been activated:

```language=bash
mamba activate moose
```

Once you're in the Gnat directory and have activated your MOOSE environment, proceed to the next step.

## 3. Getting Cardinal

Gnat allows for neutronics coupling with the continuous-energy Monte Carlo
code [OpenMC](https://github.com/openmc-dev/openmc) through
[Cardinal](https://github.com/neams-th-coe/cardinal). If you don't want to build with Cardinal,
you can skip this step. If you want to build Gnat with Cardinal, you need to update the Cardinal
submodule and fetch it's dependencies:

```language=bash
git submodule update --init cardinal
./cardinal/scripts/get-dependencies.sh
```

In addition to fetching Cardinal's dependencies, you'll need to export the following environment
variables:

```language=bash
# [REQUIRED] for builds with Cardinal
export ENABLE_CARDINAL=yes

# [REQUIRED WHEN USING THE MOOSE CONDA ENVIRONMENT] you must set the location of the
# root HDF5 directory provided by MOOSE for OpenMC to find
export HDF5_ROOT=$CONDA_PREFIX

# [REQUIRED ON SOME SYSTEMS] for some systems, libraries won't be linked properly unless
# you explicitly point this variable. We're working on a more elegant fix.
export LD_LIBRARY_PATH=$CONDA_PREFIX/lib:$LD_LIBRARY_PATH

# [OPTIONAL] it's a good idea to explicitly note that you are using MPI compiler wrappers
export CC=mpicc
export CXX=mpicxx
export FC=mpif90

# [REQUIRED WHEN RUNNING OPENMC] you will need cross section data at runtime;
# ythis variable must be set to point to a 'cross_sections.xml' file.
export OPENMC_CROSS_SECTIONS=${HOME}/cross_sections/endfb-vii.1-hdf5/cross_sections.xml
```

## 4. Building and Testing Gnat

To build Gnat, you can use the `make` buildsystem:

```language=bash
make -j{NUM_PROCESSES}
```

Where `{NUM_PROCESSES}` is the number of concurrent processes you want building
both Gnat and the MOOSE framework.
You may also compile a debug version of Gnat by running `METHOD=dbg make
-j{NUM_PROCESSES}`.

Finally, run all tests to ensure that Gnat has been built successfully:

```language=bash
./run_tests -j{NUM_THREADS}
```

Where `{NUM_THREADS}` is the number of threads you want running tests.

## 5. Run Simulations

Gnat simulations use the hierarchical input text (HIT) input deck specification
designed by INL for MOOSE-based applications. More information about HIT syntax
can be found in the [MOOSE tutorials](https://mooseframework.inl.gov/getting_started/examples_and_tutorials/tutorial01_app_development/step02_input_file.html#step-2-write-an-input-file). Gnat
input decks can be executed with the following shell command:

```language=bash
./gnat-opt -i input-file.i
```

Replace `input-file.i` with the name of your input file. Gnat can be executed in
parallel using MPI and the following shell command:

```language=bash
mpirun -n{NUM_THREADS} ./gnat-opt -i input-file.i
```

Where `{NUM_THREADS}` is the number of threads you want executing the input deck.

Now that Gnat has been successfully installed and tested, feel free to check out
the [tutorials](tutorials/index.md) to learn how to write a Gnat input deck.
