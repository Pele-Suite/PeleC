.. _EB-C9:

C9. Acoustic wave in cylindrical channel
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Case description
################

This is the case of an acoustic wave propagating in a cylindrical
channel. The geometry is a :math:`x`-direction aligned circular
channel with periodic boundary conditions. The acoustic wave
propagates along the :math:`x` direction.

The density pulse at the center of the channel is defined as:

.. math::
   \rho'(x) = \alpha \exp(-(x/\sigma)^2)

   \rho(x) = \rho_0 + \rho'(x)

   p(x) = p_0 + \rho'(x) c_s c_s

   u(x) = c_s \rho'(x) / \rho_0

The background pressure is set to :math:`p_0 = 100000 erg/cm^3`, the
background density is set to :math:`\rho_0 = 0.0014 g/cm^3`,
:math:`\alpha=10^{-6} g/cm^3`, and :math:`\sigma=10cm`. The length of
the channel is 100cm and the radius is 25cm. The simulations are
performed for :math:`t=0.000625s`. The CFL is set to 0.001 to minimize
time discretization errors.

Acoustic pulse at :math:`t=0.000625s`
#####################################

.. image:: /ebverification/C9/pulse.png
   :height: 200pt


Density profiles in the centerline at :math:`t=0.000625s`
#########################################################

.. image:: /ebverification/C9/rho.png
   :height: 300pt

:math:`L_2` error norm of density
#################################

The :math:`L_2` error norm for a quantity :math:`s` is defined as

.. math::
   \epsilon = \sqrt{ \frac{\sum_{i=1}^{n_x} (s_i^h-s_i^*)^2 }{n_x}}

where :math:`s^h` is the numerical solution, :math:`s^*` is the exact
solution, and :math:`n_x` is the number of cells in the
:math:`x`-direction.

.. image:: /ebverification/C9/error.png
   :height: 300pt

.. note::
   The second order convergence observed here is expected for this
   test case as all relevant physics happen in the direction
   perpendicular to the EB surface.


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
