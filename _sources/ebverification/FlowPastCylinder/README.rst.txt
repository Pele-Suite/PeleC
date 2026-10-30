.. _EB-FlowPastCylinder:

Flow Past a Cylinder
~~~~~~~~~~~~~~~~~~~~~~~~~~

Case description
################

This case simulates 2D flow past a cylinder for Reynolds numbers :math:`\text{Re} = 10,500`.  The geometry is defined as a rectangle with height :math:`D = 4` [cm] and length :math:`L = 12` [cm].  The inlet flow is constant in the positive x-direction with :math:`u(y) = U_{max}`, provided a Reynolds number

.. math::
   \text{Re} = \frac{\rho U_{max} D}{\mu},
   
where :math:`\rho` is the density and :math:`\mu` is the dynamic viscosity.


Domain definition
##################################
The cylinder with radius :math:`r=0.5` [cm] is centered at (0,0) on the upstream side of the domain as seen below. 

.. image:: /ebverification/FlowPastCylinder/domain.png
   :width: 450pt

Results: Re = 10
##################################

.. image:: /ebverification/FlowPastCylinder/Re10.png
   :width: 450pt

Results: Re = 500
##################################

.. image:: /ebverification/FlowPastCylinder/Re500.png
   :width: 450pt

