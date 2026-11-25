.. _EB-C10:

C10. Hagen–Poiseuille flow
~~~~~~~~~~~~~~~~~~~~~~~~~~

Case description
################

This is the standard Hagen–Poiseuille flow. The geometry is a
:math:`x`-direction aligned circular channel with imposed pressure
boundary conditions at the inflow and the outflow to enforce a driving
pressure gradient :math:`dp /dx`. For this case, the Mach number is
set to 0.1, the Reynolds number to 100, and the Prandtl number to
0.71. The background pressure is atmospheric and the background
temperature is 300K. The cylinder radius, :math:`R`, is 1cm and the
cylinder length, :math:`L`, is 12cm. The simulations are performed for
5 flowthroughs, at which point in time quantities such as kinetic
energy and momentum are converged.

The exact solution at steady state is

.. math::
   u(r) = \frac{G}{4 \mu} (R^2 - r^2)

where :math:`G = -dp/dx`, and :math:`\mu` is the dynamic viscosity.


Velocity profiles at :math:`x=6` cm (channel center)
####################################################

.. image:: /ebverification/C10/u.png
   :height: 300pt

:math:`L_2` error norm of velocity
##################################

The :math:`L_2` error norm for a quantity :math:`s` is defined as

.. math::
   \epsilon = \sqrt{ \frac{\int_{-r}^{r} (s^h-s^*)^2 \mathrm{d}r}{2 n_r}}

where :math:`s^h` is the numerical solution, :math:`s^*` is the exact
solution, and :math:`n_r` is the number of cells per radius.

.. image:: /ebverification/C10/error.png
   :height: 300pt

.. note::

   The exact solution centerline velocity is adjusted by 0.7% to
   account for compressibility effects. The observance of second order
   behavior is attributed to the fact that this is a diffusion
   dominated problem and diffusion treatment at EB surfaces is second
   order.

Time convergence of kinetic energy
##################################

.. image:: /ebverification/C10/ke_history.png
   :height: 300pt

.. note::

   This figure shows that the compressible solution reached a steady
   state. It is not expected that the total integration of the kinetic
   energy in the domain match the incompressible value for integrated
   kinetic energy :math:`K_e` because of compressibility effects.

Running study
#############

.. code-block:: bash

   paren=`pwd`
   pelec="${paren}/PeleC3d.gnu.MPI.ex"
   mpi_ranks=36

   res=( 4 8 16 32 )
   for i in "${res[@]}"
   do
       rm -rf "${i}"
       mkdir "${i}"
       cd "${i}" || exit
       cp "${paren}/inputs_3d" .
       ny="$((i*4))"
       nx="$((i*4*3))"
       srun -n ${mpi_ranks} "${pelec}" inputs_3d amr.n_cell="${nx} ${ny} ${ny}" > out
       ls -1v *plt*/Header | tee movie.visit
       cd "${paren}" || exit
   done
