# Building all Dependencies from Source

!alert! note title=tldr

All that you need to compile Gnat is:

```
cd $HOME
git clone https://github.com/OTU-Center-for-Small-Modular-Reactors/gnat.git
cd gnat
git submodule update --init
./moose/scripts/update_and_rebuild_petsc.sh
./moose/scripts/update_and_rebuild_libmesh.sh
./moose/scripts/update_and_rebuild_wasp.sh
make -j{NUM_PROCESSES}
```

If the above produces a `gnat-opt` executable, you can
jump straight to [#running]. If you were not successful with the above,
please consult the detailed instructions that follow.
!alert-end!

## 1. Cloning Gnat and MOOSE

When building Gnat and the entire dependency chain from source we recommend that you
use the MOOSE submodule included with Gnat to minimize the possibility for build errors.
Gnat and it's MOOSE submodule can be obtained with the following:

```language=bash
git clone https://github.com/OTU-Center-for-Small-Modular-Reactors/gnat.git
cd gnat
git submodule update --init
```

Once the repository has been cloned, proceed to the next step.

## 2. Building MOOSE's Dependencies

From here, you need to build all of MOOSE's dependencies. Before proceeding, please
ensure that your system meets [MOOSE's minimum requirements ](https://mooseframework.inl.gov/getting_started/installation/index.html).
If you're missing any of these packages, stop and install them first. The most common
build errors are caused by missing the python packages `packaging` and `pyaml`. Next,
you will need to set a few environment variables to ensure the build proceeds smoothly:

```bash
export CC=mpicc
export CXX=mpicxx
export FC=mpif90
export F77=mpif77
export F90=mpif90
```

From here, we need to build PETSc, libMesh, and WASP. This is done with a series of shell scripts
provided by MOOSE:

```bash
./moose/scripts/update_and_rebuild_petsc.sh
./moose/scripts/update_and_rebuild_libmesh.sh
./moose/scripts/update_and_rebuild_wasp.sh
```

If you encounter any errors at this stage, feel free to post the build errors on the Gnat discussion forums.

!alert! tip title=Building in parallel.
By default, the scripts provided by MOOSE will build in serial. To speed up the build process,
export the following environment variables:

```bash
export JOBS={NUM_PROCESSES}
export MOOSE_JOBS={NUM_PROCESSES}
export LIBMESH_JOBS={NUM_PROCESSES}
```

where `{NUM_PROCESSES}` is the number of processes you wish to use when building.

!alert-end!

## 3. Compile Gnat

At this point all of MOOSE's dependencies have been compiled, and you're ready to
build Gnat and MOOSE using the `make` buildsystem:


```language=bash
make -j{NUM_PROCESSES}
```

where `{NUM_PROCESSES}` is the number of concurrent processes you want building
both Gnat and the MOOSE framework.
You may also compile a debug version of Gnat by running `METHOD=dbg make
-j{NUM_PROCESSES}`.

Finally, run all tests to ensure that Gnat has been built successfully:

```language=bash
./run_tests -j{NUM_THREADS}
```

Where `{NUM_THREADS}` is the number of threads you want running tests.

## 4. Run Simulations id=running

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
