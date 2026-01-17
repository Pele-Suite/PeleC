#!/bin/bash

# Compile code
make TPL
make -j

# create inflow planes (with Mach 0.3)
mpirun -n 8 PeleC3d.gnu.MPI.ex run_initial_dump_planes.inp prob.M_inlet=0.3

# compile utility that combines planes into a turbulent inflow file
cd ../../../Submodules/PelePhysics/Support/TurbInflowGenerator
make -j
cd -
cp  ../../../Submodules/PelePhysics/Support/TurbInflowGenerator/PeleTurb3d.gnu.ex ./

# Assemble planes into turbulence file
./PeleTurb3d.gnu.ex type=diag_frame_planes ofile=INFLOW ifiles= $(ls -d output_initial/pltxcut*0) periodicity=0 0 0 normal=0

# Run a simulation using inflow for file
mpirun -n 8 PeleC3d.gnu.MPI.ex run_with_inflow.inp

# Create synthetic turbulence data (CSV file of cubic data)
cd  ../../../Submodules/PelePhysics/Support/TurbInflowGenerator
python gen_hit_ic.py -k0 4 -N 128
cd -
cp  ../../../Submodules/PelePhysics/Support/TurbInflowGenerator/hit_ic_4_128.dat ./

# Create turbulence file from synthetic turbulence data
./PeleTurb3d.gnu.ex type=turb_box hit_file=hit_ic_4_128.dat input_ncell=128 ofile=HIT

# Run with synthetic turbulent flutuations superimposed
mpirun -n 8 PeleC3d.gnu.MPI.ex run_with_hit.inp
