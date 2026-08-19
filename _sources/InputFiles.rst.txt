Input Files and Controls
------------------------

The input file specified on the command line is a free-format text file, one entry per row, that specifies input data processed by the AMReX ``ParmParse`` module.

This file needs to specified along with the executable as an `argv` option, for example:


::

	mpirun -np 64 ./Pele2d.xxx,yyy.ex inputs

Also, any entry that can be specified in the inputs file can also be specified on the command line; values specified on the command line override values in the inputs file, e.g.:

::

	mpirun -np 64 ./Pele2d.gnu.DEBUG.MPI.ex inputs amr.restart=sod_x_chk0030

The available options are divided into groups: those that control primarily AMReX are prefaced with `amr.` while those that are specific to Pele are prefaced with `pelec.`.

A typical input file looks something like the example below; a full list of Pele-specific input parameters are in `PeleC/Source/_cpp_parameters`.
These parameters, once read, are available in the `PeleC` object for use from c++.

::

    #------------------  INPUTS TO MAIN PROGRAM  -------------------
    # Stopping criteria: at least one must be specified, simulation
    # will stop when the first is met.

    # absolute stop time (s) for the simulation
    stop_time = 6

    # maximum number of time steps at base AMR level
    max_step = 30

    # maximum wall time (hr) after which simulation will be stopped
    max_wall_time = 1.0

    #---------------------------------------------------------------

    #------------------------
    # PROBLEM SIZE & GEOMETRY
    #-----------------------

    # flag for periodicity (here x direction is periodic)
    geometry.is_periodic = 1 0 0

    # 0 => cart, 1 => RZ  2=>spherical
    geometry.coord_sys   = 0

    # coordinates of domain's lower corner
    geometry.prob_lo     =   -0.3     0.0   0.0

    # coordinates of domain's upper corner
    geometry.prob_hi     =    0.3     0.3   0.15

    # number of cells along each direction at base level
    # (dx=dy=dz is required when an embedded boundary is used)
    amr.n_cell           =    128     64    32
    #---------------------------------------------------------------

    #---------------------------------------------------------------
    # PeleC specific inputs
    #---------------------------------------------------------------

    # >>>>>>>>>>>>>  BC KEYWORDS <<<<<<<<<<<<<<<<<<<<<<
    # Interior, UserBC, Symmetry, SlipWall, NoSlipWall
    # >>>>>>>>>>>>>  BC KEYWORDS <<<<<<<<<<<<<<<<<<<<<<

    # boundary condition at the lower face of each coordinate direction
    pelec.lo_bc       =  "Interior"  "UserBC"  "SlipWall"

    # boundary condition at the upper face of each coordinate direction
    pelec.hi_bc       =  "Interior"  "UserBC"  "SlipWall"

    #------------------------
    # TIME STEP CONTROL
    #------------------------

    pelec.cfl            = 0.5     # cfl number for hyperbolic system
    pelec.init_shrink    = 0.3     # first timestep is scaled by this factor
    pelec.change_max     = 1.1     # maximum factor by which timestep can increase
    pelec.dt_cutoff      = 5.e-20  # level 0 timestep below which we halt

    #------------------------
    # WHICH PHYSICS
    #------------------------

    pelec.do_hydro = 1               # enable hyperbolic term
    pelec.do_mol = 1                 # use method of lines (MOL)
    pelec.do_react = 0               # enable chemical reactions
    pelec.ppm_type = 2               # piecewise parabolic reconstruction type
    pelec.allow_negative_energy = 0  # flag to allow negative internal energy
    pelec.diffuse_temp = 0           # enable thermal diffusion
    pelec.diffuse_vel  = 0           # enable viscous diffusion
    pelec.diffuse_spec = 0           # enable species diffusion

    #------------------------
    # DIAGNOSTICS & VERBOSITY
    #------------------------

    # coarse time steps between computing integral of
    # conserved variables in the  domain
    # these values should stabilize at steady state
    pelec.sum_interval = 1

    pelec.v            = 1        # verbosity in PeleC cpp files
    amr.v              = 1        # verbosity in Amr.cpp
    #amr.grid_log       = grdlog  # name of grid logging file
    #---------------------------------------------------------------

    #---------------------------------------------------------------
    # AMR specific inputs
    #---------------------------------------------------------------

    #------------------------
    # REFINEMENT / REGRIDDING
    #------------------------

    amr.max_level       = 2       # maximum level number allowed
    amr.ref_ratio       = 2 2 2 2 # refinement ratio across levels
    amr.regrid_int      = 2 2 2 2 # how often to regrid
    amr.blocking_factor = 8       # block factor in grid generation
    amr.max_grid_size   = 64      # maximum number of cells per box along x,y,z

    # specify species name as flame tracer for
    #refinement purposes
    pelec.flame_trac_name = HO2


    #------------------------
    # TAGGING
    #------------------------
    tagging.denerr = 3             # density value
    tagging.dengrad = 0.01         # gradient of density value
    tagging.denratio = 1.1         # ratio of adjacent cells density
    tagging.max_denerr_lev = 3     # maximum level at which to use density for tagging
    tagging.max_dengrad_lev = 3    # maximum level at which to use density gradient for tagging
    tagging.max_denratio_lev = 3   # maximum level at which to use density ratio for tagging

    #------------------------
    # CHECKPOINT FILES
    #------------------------

    amr.checkpoint_files_output = 1
    amr.check_file              = chk    # root name of checkpoint/restart file
    amr.check_int               = 500    # number of timesteps between checkpoints

    #------------------------
    # PLOTFILES
    #------------------------

    amr.plot_files_output = 1
    amr.plot_file         = plt     # root name of plotfile
    amr.plot_int          = 100     # number of timesteps between plotfiles

    # pick which all derived variables to plot
    amr.derive_plot_vars  = pressure x_velocity y_velocity

    # we can initialize a solution from a plot file
    pelec.init_pltfile = "plt00000"

    #---------------------------------------------------------------

    #---------------------------------------------------------------
    # Embedded boundary (EB) inputs
    #---------------------------------------------------------------

    pelec.eb_isothermal = 1     # isothermal wall at EB
    pelec.eb_boundary_T = 300.  # EB wall temperature


    #------------------------
    # EB geometry
    #------------------------

    eb2.geom_type = sphere
    eb2.sphere_radius = 0.1
    eb2.sphere_center = 0.0 0.15 0.075
    eb2.sphere_has_fluid_inside = 0

    #---------------------------------------------------------------


.. note::

   It is possible to initialize a simulation using a plot file
   (e.g. `pelec.init_pltfile = "plt00000"`). It uses :math:`\rho`,
   :math:`u`, :math:`T`, :math:`Y` from a plot file to initialize a
   new state. The species in the new simulation will be taken from the
   plot file. The species that are not in the plot file will be set to
   zero. The species that are in the plot file but are not in the new
   simulation will be ignored (leading most probably to an error in
   species not summing to 1). It is therefore assumed that the
   non-zero species in the plot file used to initialize the simulation
   form a subset of the species in the simulation. The code will
   sanitize the species mass fractions to ensure that they fall within
   the right bounds. It will error out if the species are too far out
   of bounds (i.e., too far below 0, too far above 1, not summing to
   1). This check is controlled with `pelec.init_pltfile_massfrac_tol`
   and defaults to :math:`10^{-8}`. The plot file used to initialize
   must have the same domain as the new simulation, but may be coarser
   by a power of 2. If the plot file is coarser, it is required to also
   specify a value for the parameter `pelec.init_pltfile_coarse_levels`
   where the value is the number of factor of 2 coarsenings that separate
   the two grids. The init from plot file capability is based on the
   PelePhysics `PltFileManager <https://Pele-Suite.github.io/PelePhysics/Utility.html#plt-file-management>`_ utility.


Tagging criteria
~~~~~~~~~~~~~~~~

Tagging criteria are used to inform the refinement of flow features. They are added the input file using the `tagging` keyword (see the input file above). The following convention is used

- `*err`: tag cell for refinement when the value of the field exceeds this threshold value, i.e. :math:`f_{i,j,k} \geq v`, where :math:`f_{i,j,k}` is the field in cell :math:`(i,j,k)` and :math:`v` is the threshold.
- `*grad`: tag cell for refinement when the maximum difference of the field exceeds this threshold value, i.e.

.. math::
   \max(&|f_{i+1,j,k} - f_{i,j,k}|, |f_{i,j,k} - f_{i-1,j,k}|,\\
   &|f_{i,j+1,k} - f_{i,j,k}|, |f_{i,j,k} - f_{i,j-1,k}|,\\
   &|f_{i,j,k+1} - f_{i,j,k}|, |f_{i,j,k} - f_{i,j,k-1}|) \geq v

- `*ratio`: tag cell for refinement when the maximum ratio of the field (currently only supported for density) exceeds this threshold value, i.e.

.. math::
   \max(&|f_{i+1,j,k} / f_{i,j,k}|, |f_{i,j,k} / f_{i-1,j,k}|,|f_{i,j,k} / f_{i+1,j,k}|, |f_{i-1,j,k} / f_{i,j,k}|,\\
   &|f_{i,j+1,k} / f_{i,j,k}|, |f_{i,j,k} / f_{i,j-1,k}|,|f_{i,j,k} / f_{i,j+1,k}|, |f_{i,j-1,k} / f_{i,j,k}|,\\
   &|f_{i,j,k+1} / f_{i,j,k}|, |f_{i,j,k} / f_{i,j,k-1}|,|f_{i,j,k} / f_{i,j,k+1}|, |f_{i,j,k-1} / f_{i,j,k}|) \geq v

- `max_*_level`: maximum level for use of this tag (beyond this level, this tag will not be used for refinement).

The default values for tagging are defined in :code:`struct TaggingParm` in the `Tagging.H` file. Currently, the code supports tagging on density, pressure, velocity, vorticity, temperature, and volume fraction. However, additional tagging on fields can be leveraged through AMReX's tagging utility, see below for more details.

Additionally, tagging is supported for a user-specified species which can function as a "flame tracer" using the keyword `ftrac` and selecting the species with `pelec.flame_trac_name`. For example, the following text in the input file would tag cells for refinement where the HO2 mass fraction exceeded :math:`150 \times 10^{-6}` up to a maximum of 4 levels of refinement:

::

   pelec.flame_trac_name= HO2
   tagging.max_ftracerr_lev = 4
   tagging.ftracerr = 150.e-6

Users can specify their own tagging criteria in the `prob.H` of their case. An example of this is provided in the Taylor-Green regression test.

The above tagging criteria are implemented in PeleC. However, the user is encouraged to use the tagging functionality provided by AMReX and exposed in PeleC. Here are examples of how that is done:

::

   # Tag inside a box and a velocity magnitude value
   tagging.refinement_indicators = yLow magvel
   tagging.yLow.in_box_lo = -0.1  -0.52  -0.85
   tagging.yLow.in_box_hi =  3.1 -0.45    0.85

   tagging.magvel.max_level     = 2
   tagging.magvel.value_greater = 1.2e4
   tagging.magvel.field_name    = magvel

The following keys are implemented: `value_greater`, `value_less`, `vorticity_greater`, `adjacent_difference_greater`, `in_box_lo` and `in_box_hi` (to specify a refinement region), `max_level`, `start_time`, and `end_time`. The `field_name` key can be any derived or state variable.


Diagnostic Output
~~~~~~~~~~~~~~~~~

The verbosity flags `pelec.v` and `amr.v` control the extent of output related to the reacting flow solver and AMR grid printed during the simulation. When `pelec.v >= 1`, additional controls allow for fine tuning of the diagnostic output. The input flags `pelec.sum_interval` (number of coarse steps) and `pelec.sum_per` (simulation time) control how often integrals of conserved state quantities over the domain are computed and output. Additionally, if the `pelec.track_extrema` flag is set, the minima and maxima of several important derived quantities will be output whenever the integrals are output. By default, this includes the minimum and maximum across all massfractions, indicated by `massfrac`, but the `pelec.extrema_spec_name` can be set to `ALL` or an individual species name if this diagnostic for individual species is of interest.

To aid in the analysis of the diagnostic data, it can also be saved to log files. To do this, set `amr.data_log = datlog extremalog`, which will save the integrated values to `datlog` and the extrema to `extremalog`, if they are being computed based on the values of the flags described above. Additional problem-specific logs can also be created. Gridding information can also be recorded to a file specified with the `amr.grid_log` option.

Analysing the data a-posteriori can become extremely cumbersome when dealing with extreme datasets.
PeleC offers a set of diagnostics available at runtime (supported through PelePhysics) and more are under development.
Currently, the list of diagnostic contains:

* ``DiagFramePlane`` : extract a plane aligned in the 'x','y' or 'z' direction across the AMR hierarchy, writing
  a 2D plotfile compatible with Amrvis, Paraview or yt. Only available for 3D simulations.
* ``DiagPDF`` : extract the PDF of a given variable and write it to an ASCII file.
* ``DiagConditional`` : extract statistics (average and standard deviation, integral or sum) of a
  set of variables conditioned on the value of given variable and write it to an ASCII file.

When using `DiagPDF` or `DiagConditional`, it is possible to narrow down the diagnostic to a region of interest
by specifying a set of filters, defining a range of interest for a variable. Note also that for these two diagnostics,
fine-covered regions are masked. An arbitrary number of these diagnostics may be specified in a list by setting
``pelec.diagnostics`` in the input file and then specifying the diagnostic type and relevant inputs for each
diagnostic listed. See the
`PelePhysics Diagnostics documentation <https://Pele-Suite.github.io/PelePhysics/Utility.html#diagnostics>`_ for full
details on the options that must be specified for each diagnostic type.

::

   #--------------------------DIAGNOSTICS------------------------

    pelec.diagnostics = xnormP condT pdfTest

    pelec.xnormP.type = DiagFramePlane                             # Diagnostic type
    ...
