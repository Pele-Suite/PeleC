.. _EB-C3:

C3. Zero dimensional ignition with embedded boundaries
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Case description
################

This case tests the chemistry implementation in the presence of 
embedded boundaries. The geometry is the displacement volume of a 
piston-cylinder system as shown in the first figure. The gas mixture consists 
of H2, O2 and N2 at mass fractions 0.06,0.5 and 0.44, respectively at an 
initial temperature/pressure of 1500K/1 atm. The LiDryer chemical mechanism is used 
in this case. The variation of temperature 
over time that characterizes ignition delay is shown in the second figure where PeleC
solution compares well with Cantera.

The embedded boundary used in the zero dimensional ignition case
################################################################

.. image:: /ebverification/C3/setup.png
   :height: 300pt

Comparison of temperature transients between Cantera and PeleC
##############################################################

.. image:: /ebverification/C3/zerod_compare.png
   :height: 300pt
