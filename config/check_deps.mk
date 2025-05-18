define n


endef

# Set the default to the location of the Cardinal submodule, or use user-defined environment variable
CARDINAL_DIR         ?= $(CURDIR)/cardinal

# By default, don't enable Cardinal
ENABLE_CARDINAL      ?= no

# Check to ensure the user checked out the Cardinal submodule.
ifeq ($(ENABLE_CARDINAL),yes)
	CARDINAL_CONTENT     := $(shell ls $(CARDINAL_DIR) 2> /dev/null)
endif

ifeq ($(CARDINAL_CONTENT),)
  $(warning $n"Cardinal does not seem to be available. If usage of Cardinal is desired within Gnat, make sure that either the submodule is checked out$nor that CARDINAL_DIR points to a location with the Cardinal source. $n$nIn the meantime, Gnat will be built without Cardinal.")
  ENABLE_CARDINAL    := no
else
  $(info Gnat is using Cardinal from     $(CARDINAL_DIR))
endif
