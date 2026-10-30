.. highlight:: rst

.. Warning:: This documentation is a work in progress. It is reasonably complete and hopefully useful but is likely error prone and in places misleading.


Introduction
============

Pele solves the reacting compressible Navier-Stokes on a structured grid with, optionally, embedded boundary geometry treatment and non-ideal gas equations of state. A variety of time advance schemes are implemented, notably an operator-split (Strang) approach and an SDC based approach. A variety of examples are included to provide a model setup for the various options. These are discussed further in the :ref:`Getting Started<GettingStarted>` section.

The Pele Project
----------------

The origin of PeleC is in the context of *The Pele Project*, with the goal of providing a platform for combustion research in the ExaScale computing era. Specifically, to enable research grade simulation (DNS, or near DNS hybrid LES/DNS) of turbulence-chemistry interaction under conditions motivated by internal combustion engine research as well as establish a path for development of scalable design codes suitable for exascale hardware.

.. csv-table:: Design approach for Pele
   :header: "Characteristic / Need", "Approach"
   :widths: 10, 10

      "Impulsively started jets with disparate scale between fronts and turbulence; localized flames", "Dynamic adaptive mesh refinement"
      "High speed injection followed by subsonic conditions downstream", "Compressible and low-Mach capabilities"
      "Long time horizons to set up turbulence for studying fundamental TCI", "Hybrid DNS/LES [Non-reacting LES, DNS for flame]"
      "Lean, rich, and low temperature chemistry critical in multi-stage ignition and formation of soot precursors", "Accurate and detailed thermochemistry"
      "Liquid fuel injection", "Lagrangian spray model"
      "Coupling between mixture preparation and emissions", "Detailed kinetics including emissions, sectional model for soot with radiation"
      "Mixture preparation dependent on re-entrainment of combustion products", "Realistic piston dish and cylinder wall geometry"



Dependencies
------------

Pele is built on AMReX (available at `https://github.com/AMReX-Codes/amrex <https://github.com/AMReX-Codes/amrex>`_), an adaptive mesh refinement software framework, which provides the underlying software infrastructure for block structured AMR operations. The full AMReX documentation can be found `here <https://amrex-codes.github.io/amrex>`_. 

Modules for describing the equation of state, diffusion transport model, and reaction kinetics are localized to the ``PelePhysics`` repository. For the purpose of this Users' Guide  ``PelePhysics`` is considered part of ``PeleC`` but it needs to be obtained through a separate checkout as described in the :ref:`readme <README>`.


Development
-----------

A separate developers guide does not yet exist; along with the algorithmic description in this Users' Guide, doxygen documentation exists in place and an input file exists in `PeleC/Docs` that can be build using:

::

	doxygen Doxyfile
