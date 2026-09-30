# GNUmake build for the NSCBC field-dump tool, for those not using CMake.
# Modelled on amrex/Tools/Plotfile/GNUmakefile: fielddump.cpp carries its own
# main(), so the multiple_executables pattern builds it directly.
#
#   make -j              -> fielddump.<machineSuffix>.ex
#   make COMP=llvm -j    -> match whatever you built PeleC with
#
# DIM must match the plotfiles being read; fielddump.cpp is written against
# AMREX_SPACEDIM, so DIM=3 works and metrics.py detects the dimensionality from
# the header it writes.
#
#   make DIM=3 -j
#   make AMREX_HOME=/path/to/PelePhysics/Submodules/amrex DIM=3 -j

PELEC_HOME  ?= ../..
AMREX_HOME  ?= $(PELEC_HOME)/Submodules/PelePhysics/Submodules/amrex

DEBUG     = FALSE
DIM      ?= 2
COMP     ?= gnu
PRECISION = DOUBLE
USE_MPI   = FALSE
USE_OMP   = FALSE
USE_CUDA  = FALSE

BL_NO_FORT = TRUE

programs += fielddump

include $(AMREX_HOME)/Tools/GNUMake/Make.defs

multiple_executables = $(addsuffix .$(machineSuffix).ex, $(programs))
default: $(multiple_executables)

include $(AMREX_HOME)/Src/Base/Make.package
include $(AMREX_HOME)/Src/Boundary/Make.package
include $(AMREX_HOME)/Src/AmrCore/Make.package

include $(AMREX_HOME)/Tools/GNUMake/Make.rules

clean::
	$(SILENT) $(RM) $(multiple_executables)
