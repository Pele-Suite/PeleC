# Subsonic flow in a converging nozzle geometry

This simulation is of a subsonic flow in a converging nozzle
geometry. The initial condition is a quiescent flow field. The inflow
velocity is specified on the low x-direction boundary. Outflow is
specified at the high x-direction boundary using first order
extrapolation.

This example includes an optional demonstration of the "TurbInflow" capability,
which can add turbulent fluctuations to an inflow boundary or set the inflow
boundary condition based on planes saved from a previous simulation. This uses
capability from PelePhysics, see that documentation for more details on the
[TurbInflow](https://amrex-combustion.github.io/PelePhysics/Utility.html#turbulent-inflows)
capability and the [TurbInflowGenerator](https://amrex-combustion.github.io/PelePhysics/Support.html#turbfilegenerator) support code. The TurbInflow utility populates the `turb_fluc`
array that can be used to set the boundary condition in the `bcnormal` function in `prob.H`

## TurbInflow: Added Velocity fluctuations

The first step is to create the file containing the velocity fluctuation data. This is done in
two steps, first using a python code to generate a cube of synthetic homogeneous isotropic
turbulence data as a text data file and second using a C++ code to compile it into the format
expected by the PelePhysics TurbInflow utility, saved in a directory labelled `HIT`.

```bash
cd ../../../Submodules/PelePhysics/Support/TurbInflowGenerator
python gen_hit_ic.py -k0 4 -N 128
make -j
./PeleTurb3d.gnu.ex type=turb_box hit_file=hit_ic_4_128.dat input_ncell=128 ofile=HIT
cd -
cp -r ../../../Submodules/PelePhysics/Support/TurbInflowGenerator/HIT ./HIT
```

Then build and run the code using the modified inputs in `run_with_hit.inp`.

```bash
make TPL
make -j
mpirun -n 8 PeleC3d.gnu.MPI.ex turb_inflow_base.inp FILE=run_with_hit.inp
```

## TurbInflow: From saved planes from a precursor simulation

This example is contrived and is just meant to demonstrate capability. The precursor
simulation is the same geometry as the main simulation but with a different Mach number
and its inflow plane is saved and used to create the inflow for a new simulation, adjusting
the Mach number. If all the user desires is to adjust the Mach number, in practice it would
be directly changed in the inflow file. But this workflow could be used for more complex
applications like running a fully developed turbulent channel flow and using that as
an inflow for a downstream simulation.

First, run a precursor simulation with the modified inputs in `run_initial_dump_planes.inp`:

```bash
make TPL
make -j
mpirun -n 8 PeleC3d.gnu.MPI.ex turb_inflow_base.inp FILE=run_initial_dump_planes.inp
```

Then, use TurbInflowGenerator to compile the planes into the expected TurbInflow file format:

```bash
cd ../../../Submodules/PelePhysics/Support/TurbInflowGenerator
make -j
cd -
cp  ../../../Submodules/PelePhysics/Support/TurbInflowGenerator/PeleTurb3d.gnu.ex ./
./PeleTurb3d.gnu.ex type=diag_frame_planes ofile=INFLOW ifiles= $(ls -d output_initial/pltxcut*[0-9]) periodicity=0 0 0 normal=0
```

Finally, run the main simulation using the inflow generated from the precursor simulation. Note
that the main simulation duration must be fully covered by the duration of the simulation
used to generate the inflow.

```bash
mpirun -n 8 PeleC3d.gnu.MPI.ex turb_inflow_base.inp FILE=run_with_inflow.inp
```