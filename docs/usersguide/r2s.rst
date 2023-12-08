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
[2]_ and the ALARA activation code [3]_.
For Cartesian mesh, the CAD geometry must be discretized onto the mesh in order
to obtain material compositions in each mesh volume element for activation.
This is done using a ray-tracing technique. For both Cartesian and tetrahedral
meshes, mesh-based photon sampling is accomplished using MCNP5 compiled with a
custom source subroutine which utilizes the functionality of the
pyne.source_sampling module. PyNE R2S has been validated
using the Frascati Neutron Generator ITER benchmark problem, with close
agreement to experimental results [4]_.

PyNE R2S now supports DAG-OpenMC [5]_.

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
2. An MCNP MESHTAL file, DAG-MCNP tetrahedral mesh tally, or OpenMC statepoint
   file containing the neutron fluxes used for neutron activation.

Once these files have been obtained, PyNE R2S can proceed. PyNE R2S is best run
in its own folder, due the large number of intermediate files created. Create
and navigate to a new folder by executing the commands:

.. code-block:: bash

   $ mkdir my_r2s_folder
   $ cd my_r2s_folder

Then run the following command. Note that once PyNE is installed, r2s.py will
already be on $PATH.

.. code-block:: bash

   $ r2s.py setup

This command prints two new files. The first file is the configuration file 
"config.ini". Fill out this file with the appropriate information. The second
file is alara_params.txt. This file will be appended to the ALARA geometry
information automatically created by PyNE R2S to form a complete ALARA input
file. The alara_params.txt file will need to be modified to specify the irradiation schedule,
decay times of interest, and any other parameters. More information can be found in the
`ALARA Users' Guide <http://svalinn.github.io/ALARA/usersguide/index.html>`_. Once both
of these files are filled out, run the command:

.. code-block:: bash

   $ r2s.py step1

This command will generate the necessary input for running ALARA. The only
remaining files necessary will be the ALARA data file (.lib, .gam) and an
ALARA nuclib. Assuming ALARA is installed, ALARA can then be run with the 
command:

.. code-block:: bash

   $ alara alara_geom > out.txt

For large problems (i.e large meshes, many decay times), this process may take
a large amount of processor time and RAM. Once this process is complete,
execute the final command:

.. code-block:: bash

   $ r2s.py step2

This command will generate photon source density distribution meshes, one per
decay time. These files will be named like:

source_1.h5m, source_2.h5m ... source_N.h5m

A variety of information for source sampling is tagged in these *source_x.h5m*
files:

1. *e_bounds*, the energy boundaries (eV) of the photon source.
2. *decay_time*, the decay time (s) after shutdown.
3. *source_intensity*, the total photon source intensity of the specific decay
   time.
4. *r2s_source_file_version*, the version of the source file.
5. *cell_number*, the cell list in each mesh element.
6. *cell_fracs*, the cell volume fraction in the mesh element.

These source files with can now be used as sources for photon transport within
DAG-MCNP or DAG-OpenMC. Information on compiling/using a version of MCNP5 that
can utilize these mesh-based sources is found in the PyNE user's guide entry on
`mesh-based source sampling <http://pyne.io/usersguide/source_sampling.html#source-sampling-in-mcnp5>`_.
Note that each of these source files must be renamed to "source.h5m" for this
purpose. By using these sources for photon transport, the shutdown dose rate
can be obtained. Tally results will have to be normalized by the total photon
source intentity, which can be found in the "total_photon_source_intensites.txt"
file printed out by r2s.py step2 or the *source_intensity* tag in *source.h5m*.

****************
PyNE R2S example
****************

Using a simple geometry as a example, here is how we perform R2S calculation.
The example geometry is composed of four cubes of size 10x10x10 cm\ :sup:`3`\ .
It consists of two different materials, the water (in blue) and the steel
(in red). The upper left cube (in white) is void. There is an isotropic
monoenergetic (14 MeV) neutron point source in the center of the bottom right
cube. The model is irradiated by a neutron source with an intensity of :math: `1 \times 10^{10}` neutrons / second.
for 3.5 days. The following content demonstrates the process of caculating the
shutdown dose rate of this model at the time of 1 hour after shutdown.
Example files can be found in `r2s_examples <https://github.com/pyne/data/blob/master/r2s_examples.tar.gz>`_.


.. figure:: r2s_example_geometry.png
    :align: center

    **Figure 1:** *Geometry of R2S example (X-Y cross section).*

Build the model using Cubit (TM) and set material group and export to `geom_without_material.h5m`.

Prepare material library, build material library. The material library is then
generated in file "example_material_lib.h5".

.. code-block:: bash

   $ python make_example_material.py

Combine the geometry file and the material library, using  the following command:

.. code-block:: bash

   $ cp geom_without_material.h5m geom.h5m

   $ uwuw_preproc geom.h5m -v -l example_material_lib.h5

Prepare input file, define source, tally and other data cards. Example input
file can be seen in r2s_examples/neutron_transport/input.

Neutron transport calculation. A meshtal will be generated in this step.

.. code-block:: bash

   $ cp neutron_transport

   $ ln -sf ../geom.h5m .

   $ mcnp5.mpi i=input g=geom.h5m

Perform R2S setup. The 'alara_params.txt' and 'config.ini' will be generated in this step.

.. code-block:: bash

   $ cp r2s_run

   $ ln -sf ../neutron_transport/meshtal .

   $ r2s.py setup

Modify important parameters in the 'alara_params.txt' and 'config.ini' according to the problem.
There are two modes of R2S: the voxel R2S and sub-voxel R2S [6]_.
Examples input files can be seen in both 'r2s_examples/r2s_run' and 'subvoxel_r2s_run'.
By setting the parameter 'sub_voxel' (in 'config.ini') to be 'False', the user can perform voxel R2S described in [4]_.
By setting the parameter 'sub_voxel' to be 'True', the user can perform sub-voxel R2S without any other change.
Prepare alara nuclide library, copy preinstalled data library from
`ALARA/data <https://github.com/svalinn/ALARA/tree/master/data>`_.
Example nuclide library can be seen in 'r2s_examples/data'
Perform R2S step1. ALARA input file and neutron flux file will be generated in this step.

.. code-block:: bash

   $ r2s.py step1

Perform R2S step2. Several photon source files will be generated in this step.

.. code-block:: bash

   $ r2s.py step2

Perform Photon transport calculation. Example input file can be seen in r2s_examples/photon_transport/input.

.. code-block:: bash

   $ cd photon_transport

   $ ln -sf ../geom.h5m .

   $ ln -sf ../r2s_run/source_1.h5m source.h5m

   $ mcnp5.mpi i=input g=geom.h5m

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

.. [5] X. Zhang, P. C. Shriwise, S. Liu, “Implementation and Verification of
       PyNE R2S with DAG-OpenMC” Fusion Engineering and Design, Vol. 160:
       111837, (2020).

.. [6] X. ZHANG, P. C. Shriwise, S. Liu, P. P.H. Wilson, “PyNE Sub-Voxel R2S
       for Shutdown Dose Rate Analysis”, International Conference on
       Mathematics, Computational Methods & Reactor Physics (M&C 2019), (2019).
