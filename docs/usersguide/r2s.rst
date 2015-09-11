.. _usersguide_r2s:

==================================
Rigorous Two-Step Activation (R2S)
==================================

.. currentmodule:: pyne.r2s

.. automodule:: pyne.r2s

********
Overview
********

The Rigorous Two-step method [1] is a method for estimating the photon dose
rates that result from neutron activation, often as a function of position and
time after "shutdown" (i.e. when neutron production ceases). This so-called
"shutdown dose rate is calculated using two seperate transport steps, using the procedure below:

1. Neutron transport to obtain a global, energy-wise neutron flux distribution, typically on a mesh.
2. Nuclear inventory analysis to calculate to an energy photon emission density
   distribution for each time after shutdown of interest.
3. Photon transport using each of the photon emission densities found in 2. as sources in order to
   calculate photon dose rates.

PyNE R2S impliments Carteian- and tetrahedral- mesh-based R2S method that
operated entirely on CAD mesh. For Cartesian mesh, the geometry. This is
accomplished using the using the Directed Accelerated Geometry Monte Carlo
(Version) of MCNP5, known as DAGMC5 and the ALARA activation code.



**********************
Command-Line Interface
**********************

The source sampling module implements mesh-based source sampling

located in the scripts/ directory. 

.. code-block:: bash

   >> r2s.py setup

.. code-block:: bash

   >> r2s.py step1


.. code-block:: bash

   >> alara alara_geom > out.txt


.. code-block:: bash

   >> r2s.py step2

.. code-block:: bash

  idum 1 100

:analog:
  Particle birth parameters are sampled directly from a unmodified 


**********
References
**********

[1] Y. Chen, U. Fischer, Rigorous MCNP Based Shutdown Dose Rate Calculations:
Computational Scheme, Verification Calculations and Application to ITER, Fusion
Engineering and Design, Vol. 63-64, (2002)

[2] E. Biondo, A. Davis, A. Scopatz, P. Wilson, "Rigorous Two-Step Activation for
Fusion Systems with PyNE", Transactions of the American Nuclear Society,
Vol. 112, (2015).

Mesh-Based Source Sampling http://pyne.io/usersguide/source_sampling.html

[3] The University of Wisconsin Unified Workflow, http://svalinn.github.io/DAGMC/usersguide/uw2.html

[4] ALARA Users' Guide. http://alara.engr.wisc.edu/users.guide.html/

