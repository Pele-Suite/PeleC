.. _EB-C7:

C7. Sod shock tube in rotated channel with AMR
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Case description
################

This is the standard Sod shock tube problem. The geometry is a
rectangular channel at :math:`30^\circ` to the mesh. This case is run
with multiple levels of AMR. The profiles shown below are through the
channel centerline.

Density at t=0.1s
#################

.. image:: /ebverification/C7/stube.png
   :height: 300pt

Field profiles in the centerline at t=0.1s
##########################################

.. image:: /ebverification/C7/rho.png
   :height: 300pt

.. image:: /ebverification/C7/p.png
   :height: 300pt

.. image:: /ebverification/C7/u.png
   :height: 300pt

.. image:: /ebverification/C7/temp.png
   :height: 300pt


Running study
#############

.. code-block:: bash

   paren=`pwd`
   pelec="${paren}/PeleC3d.gnu.MPI.ex"
   mpi_ranks=36

   lev=( 0 1 2 )
   for i in "${lev[@]}"
   do
       rm -rf "${i}"
       mkdir "${i}"
       cd "${i}" || exit
       cp "${paren}/inputs_3d" .
       srun -n ${mpi_ranks} "${pelec}" inputs_3d `python3 ${paren}/gen_tube_input.py -a 30.0 -n 8` amr.max_level="${i}" > out
       ls -1v *plt*/Header | tee movie.visit
       cd "${paren}" || exit
   done
