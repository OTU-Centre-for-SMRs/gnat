# These modules are necessary for Gnat.
NAVIER_STOKES               := yes
RAY_TRACING                 := yes

# If building with Cardinal, the Reactor module is necessary.
ifeq ($(ENABLE_CARDINAL),yes)
	REACTOR                     := yes
else
	REACTOR                     ?= yes
endif
