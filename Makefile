# Whether or not Gnat should be built with Cardinal.
ENABLE_CARDINAL      ?= no

###############################################################################
################### MOOSE Application Standard Makefile #######################
###############################################################################
#
# Optional Environment variables
# MOOSE_DIR        - Root directory of the MOOSE project
#
###############################################################################
# Use the MOOSE submodule if it exists and MOOSE_DIR is not set
MOOSE_SUBMODULE    := $(CURDIR)/moose
ifneq ($(wildcard $(MOOSE_SUBMODULE)/framework/Makefile),)
  MOOSE_DIR        ?= $(MOOSE_SUBMODULE)
else
  MOOSE_DIR        ?= $(shell dirname `pwd`)/moose
endif
FRAMEWORK_DIR      := $(MOOSE_DIR)/framework

include config/check_deps.mk

ifeq ($(ENABLE_CARDINAL),yes)
  include config/configure_cardinal.mk
endif

include $(FRAMEWORK_DIR)/build.mk
include $(FRAMEWORK_DIR)/moose.mk

################################## MODULES ####################################
# Set the mandatory modules.
include config/set_modules.mk

# These modules are optional.
FLUID_PROPERTIES            := yes
HEAT_TRANSFER               := yes
THERMAL_HYDRAULICS          := yes

include $(MOOSE_DIR)/modules/modules.mk
###############################################################################

# Build Cardinal
ifeq ($(ENABLE_CARDINAL),yes)
	include config/build_cardinal.mk
endif

# dep apps
APPLICATION_DIR    := $(CURDIR)
APPLICATION_NAME   := gnat
BUILD_EXEC         := yes
GEN_REVISION       := yes

# Cardinal dependency libraries needed for linking
ifeq ($(ENABLE_CARDINAL),yes)
  ADDITIONAL_LIBS := -L$(CARDINAL_DIR)/lib $(CC_LINKER_SLFLAG)$(CARDINAL_DIR)/lib \
                     -L$(OPENMC_LIBDIR) -lopenmc -lhdf5_hl -ldagmc -lMOAB \
                     $(CC_LINKER_SLFLAG)$(OPENMC_LIBDIR)
endif

include            $(FRAMEWORK_DIR)/app.mk

###############################################################################
# Additional special case targets should be added here
ifeq ($(ENABLE_CARDINAL),yes)
  # Cardinal contrib flags used in app.mk targets
  CARDINAL_EXTERNAL_FLAGS := -L$(CARDINAL_DIR)/lib $(CC_LINKER_SLFLAG)$(CARDINAL_DIR)/lib $(BLASLAPACK_LIB) $(PETSC_EXTERNAL_LIB_BASIC) -L$(OPENMC_LIBDIR)  -L$(HDF5_LIBDIR) -lopenmc -ldagmc -lMOAB $(CC_LINKER_SLFLAG)$(OPENMC_LIBDIR) $(CC_LINKER_SLFLAG)$(HDF5_LIBDIR)

  # EXTERNAL_FLAGS is used in rules for app.mk
  $(app_LIB): EXTERNAL_FLAGS := $(CARDINAL_EXTERNAL_FLAGS)
  $(app_test_LIB): EXTERNAL_FLAGS := $(CARDINAL_EXTERNAL_FLAGS)
  $(app_EXEC): EXTERNAL_FLAGS := $(CARDINAL_EXTERNAL_FLAGS)
endif
