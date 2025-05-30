# Add OpenMC flags
ADDITIONAL_CPPFLAGS += $(HDF5_INCLUDES) $(OPENMC_INCLUDES)
libmesh_CXXFLAGS    += -DENABLE_OPENMC_COUPLING

# Add DagMC flags (-DDAGMC is used in OpenMC)
libmesh_CXXFLAGS    += -DENABLE_DAGMC -DDAGMC

# For some reason we need to define these...
# TODO: fix this.
libmesh_CXXFLAGS    += -DCARDINAL_INSTALLABLE_DIRS=GNAT_INSTALLABLE_DIRS
libmesh_CXXFLAGS    += -DCARDINAL_REVISION=GNAT_REVISION
libmesh_CXXFLAGS    += -DCARDINAL_VERSION=GNAT_VERSION

# Configure and build MOAB, DagMC, and then OpenMC
include          $(CARDINAL_DIR)/config/moab.mk
include          $(CARDINAL_DIR)/config/dagmc.mk

# autoconf-archive puts some arguments (e.g. -std=c++17) into the compiler
# variable rather than the compiler flags variable.
#
# cmake allows this, but wants any compiler arguments to be
# semicolon-separated, not space-separated
# libmesh_CC, etc., were defined in build.mk
space := $(subst ,, )
LIBMESH_CC_LIST := $(subst $(space),;,$(libmesh_CC))
LIBMESH_CXX_LIST := $(subst $(space),;,$(libmesh_CXX))
LIBMESH_F90_LIST := $(subst $(space),;,$(libmesh_F90))

ENABLE_DAGMC     := ON
include            $(CARDINAL_DIR)/config/openmc.mk

# Cardinal
libmesh_CXXFLAGS   += -DENABLE_CARDINAL
APPLICATION_DIR    := $(CARDINAL_DIR)
APPLICATION_NAME   := cardinal
BUILD_EXEC         := no
GEN_REVISION       := yes
include            $(FRAMEWORK_DIR)/app.mk

# app_objects are defined in moose.mk and built according to the rules in build.mk
# We need to build these first so we get include dirs
$(app_objects): build_moab build_dagmc build_openmc
$(test_objects): build_moab build_dagmc build_openmc
