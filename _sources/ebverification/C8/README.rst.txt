.. _EB-C8:

C8. Multi-species shock tube in a rotated channel with AMR
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Case description
################

This case is a non-reacting multi-species shock tube problem 
that tests inviscid hydrodynamics and the handling 
of species in a domain with embedded boundaries. The geometry is a
rectangular channel at :math:`30^\circ` to the mesh. The left half of the 
channel is initialized with pure nitrogen while the right is initialized with helium.
The pressure and density ratios are 0.1 and 0.125 respectively, similar to the Sod shock tube problem.
This case is run with multiple levels of AMR. The figures shown below indicate mesh refinement at 
the embedded boundary, contact and shock discontinuities. The center-line density profile is 
compares well with data from literature.

AMR mesh indicating refinement at embedded boundary, contact and shock discontinuities
######################################################################################

.. image:: /ebverification/C8/ebmesh.png
   :height: 300pt

Field profiles in the centerline at t=0.2 L sqrt(rhoL/pL)
#########################################################

.. image:: /ebverification/C8/density.png
   :height: 300pt
