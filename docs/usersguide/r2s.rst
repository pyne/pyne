.. _usersguide_r2s:

==================================
Rigorous Two-Step Activation (R2S)
==================================

.. currentmodule:: pyne.r2s

.. automodule:: pyne.r2s

********
Overview
********

The Rigorous Two-Step (R2S) method [1]_ is a method for estimating the photon dose
rates that result from neutron activation, often as a function of position and
time after "shutdown" (i.e. when neutron production ceases). The primary
application of this method is occupational safety and maintenance planning
within facilities that generate neutrons such as fission reactors, fusion
devices, and experimental facilities. The so-called "shutdown dose rate" is
calculated using two separate transport steps, using the procedure below:

1. Neutron transport to obtain a global, energy-wise neutron flux distribution,
   typically on a mesh.
2. Nuclear inventory analysis to calculate an energy-wise photon emission density
   distribution for each time after shutdown of interest.
3. Photon transport using each of the photon emission density distributions found as
   sources in order to calculate photon dose rates.

PyNE R2S implements a Cartesian- and tetrahedral- mesh-based R2S method that
operates entirely on CAD geometry. This is accomplished using the Direct
Accelerated Geometry Monte Carlo (DAGMC) version of MCNP5, known as DAG-MCNP5
[2]_ and the ALARA activation code [3]_. For Cartesian mesh, the CAD geometry
must be discretized onto the mesh in order to obtain material compositions in each
mesh volume element for activation. This is done using a ray-tracing technique.
For both Cartesian and tetrahedral meshes, mesh-based photon sampling is
accomplished using MCNP5 compiled with a custom source subroutine which utilizes the
functionality of the pyne.source_sampling module. PyNE R2S has been validated
using the Frascati Neutron Generator ITER benchmark problem, with close
agreement to experimental results [4]_.

***************
Using PyNE R2S
***************

PyNE R2S is principally used through the command-line interface (CLI) found
in pyne/scripts/r2s.py. This section explains how to use PyNE R2S via this
script. Note that more complex use cases may not be handled by this script,
in which case the user can use the python interface in pyne/pyne/r2s.py. In
order to use the CLI, the following files are first needed:

1. A material-laden CAD file representing the geometry. Instructions on creating
   this file can be found `here <http://svalinn.github.io/DAGMC/usersguide/uw2.html>`_.
2. An MCNP MESHTAL file or DAG-MCNP tetrahedral mesh tally containing the
   neutron fluxes used for neutron activation.

Once these files have been obtained, PyNE R2S can proceed. PyNE R2S is best run
in its own folder, due the large number of intermediate files created. Create
and navigate to a new folder by executing the commands:

.. code-block:: bash

   >> mkdir my_r2s_folder
   >> cd my_r2s_folder

Then run the following command. Note that once PyNE is installed, r2s.py will
already be on $PATH.

.. code-block:: bash

   >> r2s.py setup

This command prints two new files. The first file is the configuration file 
"config.ini". Fill out this file with the appropriate information. The second
file is alara_params.txt. This file will be appended to the ALARA geometry
information automatically created by PyNE R2S to form a complete ALARA input
file. The alara_params.txt file will need to be modified to specify the irradiation schedule,
decay times of interest, and any other parameters. More information can be found in the
`ALARA Users' Guide <http://alara.engr.wisc.edu/users.guide.html/>`_. Once both
of these files are filled out, run the command:

.. code-block:: bash

   >> r2s.py step1

This command will generate the necessary input for running ALARA. The only
remaining files necessary will be the ALARA data file (.lib, .gam) and an
ALARA nuclib. Assuming ALARA is installed, ALARA can then be run with the 
command:

.. code-block:: bash

   >> alara alara_geom > out.txt

For large problems (i.e large meshes, many decay times), this process may take
a large amount of processor time and RAM. Once this process is complete,
execute the final command:

.. code-block:: bash

   >> r2s.py step2

This command will generate photon source density distribution meshes, one per
decay time. These files will be named like:

source_1.h5m, source_2.h5m ... source_N.h5m

These files can now be used as sources
for photon transport within MCNP. Information on compiling a version of MCNP5
that can utilize these sources is found in the PyNE user's guide entry on 
`mesh-based source sampling <http://pyne.io/usersguide/source_sampling.html#source-sampling-in-mcnp5>`_.
Note the each of these source files must be renamed to "source.h5m" for this purpose. By using these sources for photon transport, the shutdown dose rate can be obtained as a final answer.

**********
References
**********

.. [1] Y. Chen, U. Fischer, Rigorous MCNP Based Shutdown Dose Rate Calculations:
       Computational Scheme, Verification Calculations and Application to ITER, Fusion
       Engineering and Design, Vol. 63-64, (2002)

.. [2] T. J. Tautges, P.P.H.  Wilson, J. Kraftcheck, B. F. Smith, D. L. 
       Henderson, "Acceleration Techniques for Direct Use of CAD-Based 
       Geometries in Monte Carlo Radiation Transport", International Conference
       on Mathematics, Computational Methods & Reactor Physics (M&C 2009),
       (2009).

.. [3]  P. P.H. Wilson, H. Tsige-Tamirat, H. Y. Khater, D. L.  Henderson,
        Validation of the ALARA Activation Code, Fusion Technology, Vol. 34,
        Issue 3, (1998),

.. [4] E. Biondo, A. Davis, A. Scopatz, P. P.H. Wilson, "Rigorous Two-Step 
       Activation for Fusion Systems with PyNE", Transactions of the American 
       Nuclear Society, Vol. 112, (2015).

