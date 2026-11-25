.. _EB-C1:

C1. Method of manufactured solutions for EB
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Case description
################

This is a similar test case problem as performed for PeleC
verification using MASA to construct a manufactured solution. The MMS
is performed around a sphere at a Re number and Mach number of 1.

:math:`L_2` error norm of fields
################################

.. image:: /ebverification/C1/rho_error.png
   :height: 300pt

.. image:: /ebverification/C1/p_error.png
   :height: 300pt

.. image:: /ebverification/C1/u_error.png
   :height: 300pt

.. image:: /ebverification/C1/v_error.png
   :height: 300pt

.. image:: /ebverification/C1/w_error.png
   :height: 300pt

.. note::
   The first order convergence observed here is expected because the
   treatment of the EB surface is first order.


Running study
#############

.. code-block:: bash

   paren=`pwd`
   pelec="${paren}/PeleC3d.gnu.MPI.ex"
   mpi_ranks=36

   res=( 8 16 32 64 )
   for i in "${res[@]}"
   do
       rm -rf "${i}"
       mkdir "${i}"
       cd "${i}" || exit
       cp "${paren}/inputs_3d" .
       srun -n ${mpi_ranks} "${pelec}" inputs_3d amr.n_cell="${i} ${i} ${i}" > out
       ls -1v *plt*/Header | tee movie.visit
       cd "${paren}" || exit
   done
