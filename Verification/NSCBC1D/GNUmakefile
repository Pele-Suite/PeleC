# GNUmake build for the standalone 1-D NSCBC driver, for those not using CMake.
#
# The CMakeLists.txt alongside this file needs an INSTALLED AMReX
# (-DAMReX_DIR=<prefix>/lib/cmake/AMReX).  If what you have is a source tree --
# a PelePhysics submodule checkout, say -- use this instead:
#
#   make -j                                  # air, 2 species: 22 checks
#   make Chemistry_Model=LiDryer -j          # a reacting mechanism: 26 checks,
#                                            # C7 (the reaction source) included
#
# AMREX_HOME may come from the environment, which is the usual thing when the
# amrex you build against lives under another Pele repository:
#
#   export AMREX_HOME=$HOME/src/PeleLMeX/Submodules/PelePhysics/Submodules/amrex
#   make -j && ./nscbc1d.gnu.ex
#
# `./nscbc1d.gnu.ex sweep` also dumps the sigma table.
#
# Like the CMake build this compiles exactly two translation units: the driver
# and the mechanism.  PelePhysics' EOS is header-only at this level, so none of
# its .cpp files are needed and none of its heavier dependencies are pulled in.
#
# This compiles Source/NSCBC.H UNMODIFIED -- that is the whole point of the
# driver, so do not be tempted to keep a local copy of the kernel here.

PELE_HOME         ?= ../..
PELE_PHYSICS_HOME ?= $(PELE_HOME)/Submodules/PelePhysics
AMREX_HOME        ?= $(PELE_PHYSICS_HOME)/Submodules/amrex

DEBUG     = FALSE
DIM      ?= 1
COMP     ?= gnu
PRECISION = DOUBLE
USE_MPI   = FALSE
USE_OMP   = FALSE
USE_CUDA  = FALSE

BL_NO_FORT = TRUE

Chemistry_Model ?= air
PELE_EOS        ?= Fuego

programs += nscbc1d

include $(AMREX_HOME)/Tools/GNUMake/Make.defs

multiple_executables = $(addsuffix .$(machineSuffix).ex, $(programs))
default: $(multiple_executables)

include $(AMREX_HOME)/Src/Base/Make.package

CEXE_sources += mechanism.cpp

PP := $(PELE_PHYSICS_HOME)
PP_DIRS := $(PP)/Source \
           $(PP)/Source/Eos \
           $(PP)/Source/Transport \
           $(PP)/Source/Utility \
           $(PP)/Source/Utility/Utilities \
           $(PP)/Source/Utility/Filter \
           $(PP)/Source/Utility/BlackBoxFunction \
           $(PP)/Source/Utility/PMF \
           $(PP)/Source/Utility/TurbInflow \
           $(PP)/Source/Utility/Diagnostics \
           $(PP)/Mechanisms/$(Chemistry_Model)

# PeleC's Source/ carries NSCBC.H, IndexDefines.H and Constants.H.
VPATH_LOCATIONS   += . $(PELE_HOME)/Source $(PP_DIRS)
INCLUDE_LOCATIONS += . $(PELE_HOME)/Source $(PP_DIRS)

# The same set the CMake build defines: no advected or auxiliary variables.
# SRK forbids constant transport (PelePhysicsConstraints.H); pair it with
# simple transport, as the CMake build does.
ifeq ($(PELE_EOS),SRK)
DEFINES += -DUSE_SRK_EOS -DUSE_SIMPLE_TRANSPORT
else
DEFINES += -DUSE_FUEGO_EOS -DUSE_CONSTANT_TRANSPORT
endif
DEFINES += -DNUM_ADV=0 -DNUM_AUX=0

include $(AMREX_HOME)/Tools/GNUMake/Make.rules

clean::
	$(SILENT) $(RM) $(multiple_executables)
