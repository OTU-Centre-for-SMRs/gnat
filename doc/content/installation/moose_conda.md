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
git submodule update --init
```

For users who want to develop/use other MOOSE applications:

```language=bash
git clone https://github.com/OTU-Center-for-SMRs/gnat.git
```

Once the repository has been cloned, proceed to the next step.

## 3. Compile Gnat

Before continuing with the compilation process, ensure that the MOOSE environment
has been activated in your shell. This can be done with the following command:

```language=bash
mamba activate moose
```

You can then compile the application using the `make` buildsystem:

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

## 4. Run Simulations

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
