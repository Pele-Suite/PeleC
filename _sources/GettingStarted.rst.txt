 .. role:: cpp(code)
    :language: c++
 
.. _GettingStarted:

Getting Started
===============

Navigation
----------

The PeleC directory structure is as shown below:

* **Source** - C++ source code

* **Docs**   - PeleC documentation

* **Exec** - `regression tests<VandV>` and various capability demonstrations

Setting up a problem to run with PeleC involves writing an input file and problem specific code in the run directory. 
PeleC is built using the AMReX build system which supports out-of-source builds but as configured in Pele requires a specific directory structure. 
Within each case directory in Exec, are the source files that specify the setup of that particular case. 
The user has to build each case by compiling source files using a GNUMakefile which also compiles and links together AMReX and PeleC sources.
The source files contained in the case directory are treated preferentially and can override PeleC/AMReX source files. A few key files that need to be supplied for (most) cases are:

* **inputs,inp** -- a text file containing parameters that are ready by the ParmParse capability in AMReX. These include things like number of time steps, grid size, output file frequency, which physics to include, etc. A list of available data in the Pele group can be found in PeleC/Source/param_includes/pelec_params.H.

* **prob.cpp** -- Routines called at:

  * Initialization (``amrex_probinit``) 
  * Problem teardown (``pc_prob_close``)

* **prob.H** -- prob.H header file is normally used to define the user-defined embedded boundary class, if any. It is also used to define the solution initialisation and boundary condition implementation functions. User-defined functions defined in this header file include:

  * ``pc_initdata``: GPU kernel that sets initial conditions at every cell individually. Must always be defined by user.
  * ``bcnormal``: GPU kernel that sets values in ghost cells at domain boundaries. At least a dummy version must be provided. This function is called for any boundaries specified as ``UserBC`` or ``Hard``.
  * ``ProblemSpecificFunctions``: a set of problem specific capabilities that can optionally be specified by the user. To simply use the defaults,
    point this to the default class defined in ``PeleC/Source/ProblemSpecificFunctions.H`` with
    ``using ProblemSpecificFunctions = DefaultProblemSpecificFunctions;``, otherwise create a new class that inherits from the default class and
    provides implementations for the desired functions (see the TG RegTest for example). The functions and call signatures that can be defined
    are found in ``PeleC/Source/ProblemSpecificFunctions.H`` and summarized below:
    * ``add_problem_derive``: Create custom derived variables
    * ``set_problem_tags``: Create custom tagging criteria
    * ``set_aux_names``: Name the auxiliary variables if being used
    * ``set_adv_names``: Name the advected variables if being used
    * ``set_isothermal_wall_temp``: Set a spacially varying temperature for domain boundary walls (if relevant and desired)
    * ``problem_modify_ext_sources``: Add custom source terms
    * ``problem_modify_transport_coeffs``: Modify the transport coefficients (viscosity, thermal conductivity, diffusivity)

* **prob_parm.H** -- Used to define problem specific parameters and their values.

* **GNUMakefile** -- In addition to setting options to build profiling, debugging, MPI, OpenMP, Compiler toolchain options, the chemical mechanism, transport model, equation of state model, and use of EB are set here for compile time selection. The GNUMakefile includes the ``Make.PeleC`` file from the `Exec` directory that contains build configuration common across the examples.


.. toctree::
   :maxdepth: 1

   Building
   InputFiles
   Tutorials
   Testing
