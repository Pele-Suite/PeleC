
 .. role:: cpp(code)
    :language: c++

 .. _LES:

LES and Hybrid LES/DNS Support
------------------------------

.. note:: PeleC has support for LES subgrid models for *non-reacting* flows; documentation is a work in progress.


Models for large eddy simulations
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. warning:: The LES source terms do not currently support EB cut cells.


PeleC currently supports two LES models, the constant and dynamic
Smagorinsky models. An extensive discussion of the compressible
version of these models can be found in Martín, M. Pino, U. Piomelli,
and G. V. Candler. "Subgrid-Scale Models for Compressible Large-Eddy
Simulations." Theoretical and Computational Fluid Dynamics 13, no. 5
(2000): 361–76. The constant Smagorinsky model was verified using the
method of manufactured solutions.


Using
#####

To use the LES models, the user must first activate the LES source
terms by doing: ``pelec.do_les = 1``. The default LES model is
constant Smagorinsky. The user can pick the LES model by setting
``pelec.les_model = NUM``:

* ``les_model = 0``: constant Smagorinsky model
* ``les_model = 1``: dynamic Smagorinsky model
* ``les_model = 2``: WALE model
* ``les_model = 3``: Vreman model

For the constant Smagorinsky model, the user may define the model
coefficients: ``pelec.Cs``, ``pelec.CI``, and ``pelec.PrT``. These
must be set by the user in the input file as they default to ``Cs =
0``, ``CI = 0``, and ``PrT = 1.0``. To use a typical constant
Smagorinsky model, the user would add the following in the input file:

::

   pelec.do_les = 1
   pelec.les_model = 0
   # TURBULENCE PARAMETERS
   pelec.Cs = 0.16
   pelec.CI = 0.09
   pelec.PrT = 0.7


For the dynamic Smagorinsky model, the user can define a test filter
type (``les_test_filter_type``) and filter-grid ratio
(``les_test_filter_fgr``). A complete description of the filter types
and filter-grid ratios can be found in the section on filtering
below. The default test filter type is a box filter with a filter-grid
ratio of 2. An example for activating the dynamic Smagorinsky model is:

::

   pelec.do_les = 1
   pelec.les_model = 1
   pelec.les_test_filter_type = 3
   pelec.les_test_filter_fgr = 2


Explicit filtering of the hydrodynamic source terms
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. warning:: Explicit filtering is a work in progress. It should not be used with wall boundary conditions. Filtering at coarse-fine boundaries is experimental.

In some instances the user may want to explicitly filter the
Navier-Stokes equations to perform explicitly filtered large eddy
simulations. The objective of explicitly filtering the PDEs is to
remove wavenumber content in the solution that may be contaminated by
numerical error. To achieve this objective, the hydrodynamic source
terms, i.e. the source terms that may introduce wavenumbers that
cannot be represented accurately by the grid, are explicitly filtered
after each computation of the hydrodynamic source terms.

Using
#####

To explicitly filter the hydrodynamic source terms, the following
should option should be turned on in the input file:
``pelec.use_explicit_filter = 1``. The user specifies the filter-grid
ratio using ``pelec.les_filter_fgr = NUM``, where ``NUM`` is the
filter-grid ratio desired, e.g. ``pelec.les_filter_fgr = 2``. The user
also specifies a filter type through ``pelec.les_filter_type = NUM``.
Explicit filtering utilizes the Filter utility from PelePhysics, and
the set of available filter types may be found in the
`PelePhysics documentation <https://Pele-Suite.github.io/PelePhysics/Utility.html#filter>`_.

An example input file section for a Gaussian filter with a filter-grid
ration of 2 would be:

::

   pelec.use_explicit_filter=1
   pelec.les_filter_type=2
   pelec.les_filter_fgr=2
