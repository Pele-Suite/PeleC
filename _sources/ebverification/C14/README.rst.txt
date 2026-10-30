.. _EB-C14:

C14. Shock interacting with a ramp
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Case description
################

This is a shock that moves up an inclined plane. The domain size is
:math:`[0,0,0] \times [7.5, 3.75, 0.1171875]` with :math:`512 \times
256 \times 8` cells and 1 level of AMR.


Pressure
########

.. image:: /ebverification/C14/p.png
   :height: 300pt

Running study
#############

.. code-block:: bash

   paren=`pwd`
   pelec="${paren}/PeleC3d.gnu.MPI.ex"
   mpi_ranks=36

   srun -n ${mpi_ranks} "${pelec}" inputs_3d > out
   ls -1v *plt*/Header | tee movie.visit
