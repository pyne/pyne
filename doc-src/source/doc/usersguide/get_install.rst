Getting and Installing the DAGMC Toolkit
----------------------------------------

DAGMC is a toolkit that provides direct geometry support to any Monte
Carlo radiation transport code.  The primary code used for development
and testing of DAGMC is `MCNP5 <http://mcnp-green.lanl.gov/>`_,
developed by `Los Alamos National Laboratory <http://www.lanl.gov>`_
and distributed by the `Radiation Safety Information Computing Center
<http://rsicc.ornl.gov>`_.  There has also been experience with MCNPX
(LANL), Tripoli4 (CEA/Saclay), and `FLUKA <http://www.fluka.org/>`_ (CERN/INFN).

These instructions describe the basic steps for downloading and
installing the software libraries that provide the DAGMC toolkit for
integration with Monte Carlo codes.  After this, code-specific
instructions will be give for each code. 

Toolkit Dependencies and Workflows
+++++++++++++++++++++++++++++++++++++++

It is useful to consider how users will use the DAGMC workflow prior
to making some installation decisions.  There are three main stages
for the workflow:
* manual preparation of geometric models
* automated pre-processing of models
* use of the models in analysis

While the second and third stages can be combined, these instructions
will be based on treating them as separate stages with comments about
how they would be combined where appropriate.

Geometry
~~~~~~~~~~
*Manual Geometry Preparation*

Because this stage is highly interactive, most users prefer to perform
this stage on their desktop or laptop computer.  The only software
required for this stage is the interactive Cubit software.

Pre-Processing
~~~~~~~~~~~~~~~
*Automated Model Pre-processing*

The primary purpose of this stage is to "translate" from the solid
model representation to the faceted representation.  As such, it
requires Cubit, CGM and MOAB, and typically results in dependencies on
shared libraries.

Analysis
~~~~~~~~~
*Model Analysis*

During analysis, only the faceted model is necessary, requiring only
MOAB and the physics code.  It is also possible to combine this stage
with the previous one, requiring Cubit, CGM, MOAB and the physics
code, and will also result in dependencies on shared libraries.

Summary
~~~~~~~
*Summary of Dependencies* 

* Some physics package, e.g. MCNP5
   * `MOAB/DAGMC <http://trac.mcs.anl.gov/projects/ITAPS/wiki/MOAB>`_
       * `HDF5 <http://www.hdfgroup.org/HDF5/release/obtain5.html>`_
       * `CGM <http://trac.mcs.anl.gov/projects/ITAPS/wiki/CGM>`_ 
            * ACIS v19, or `CUBIT <http://cubit.sandia.gov>`_ v12.2 or v13.1 (with CGM `trunk <http://ftp.mcs.anl.gov/pub/fathom/cgm-nightly-trunk.tar.gz>`_ only)



Toolkit Installation
++++++++++++++++++++++++++++

Installing the MOAB library can be accomplished in either of two ways:

* Compile each component individually in a 4-step process.
* Run the *build_dagmc_stack.bash* script.

**Compile each component individually**

The following 4 steps are required to install the MOAB library,
including the DAGMC toolkit, for use in Monte Carlo radiation
transport tools.  Following these steps, all the pre-requisites for
the automated processing stage will be available.

1. Install `CUBIT <http://cubit.sandia.gov>`_ v12.2 or v13.1
2. Install `CGM <http://trac.mcs.anl.gov/projects/ITAPS/wiki/CGM>`_, using the --with-cubit options
3. Install `HDF5 <http://www.hdfgroup.org/HDF5/>`_
4. Install `MOAB <http://trac.mcs.anl.gov/projects/ITAPS/wiki/MOAB>`_,
   using the --with-cgm and --with-hdf5 options (--with-netcdf may
   also be useful but not necessary)

Here are some assumptions/conventions that are used in these instructions:

* all operations are carried out in the a subdirectory ``dagmc_bld`` of a user's home directory
* path to CUBIT files is known, e.g. ``/path/to/cubit``.  This is the directory that contains the Python script file ``cubit`` and a ``bin`` subdirectory.  
* all tarballs reside in user's home directory

If these do not apply to you, please modify your steps accordingly.

     *(For a shortcut to installing DAG-MCNP5.1.51 you may be able to use the DagmcBuildPackage .)*

CGM
````

CGM:  Create a directory to build in:
::
    prompt%> mkdir -p $HOME/dagmc_bld/CGM/bld
    prompt%> cd $HOME/dagmc_bld/CGM

If installing from SVN repository (you *must* use this method if you are using Cubit v13.1):
::
    prompt%> svn co https://svn.mcs.anl.gov/repos/ITAPS/cgm/trunk
    prompt%> cd trunk
    prompt%> autoreconf -fi
    prompt%> cd ..
    prompt%> ln -s trunk src

If installing CGM version 12.2 from a tarball, ``CGM-12.2.0.tar.gz``:
::
    prompt%> tar xzf ~/CGM-12.2.0.tar.gz
    prompt%> ln -s CGM-12.2.0 src

In all CGM cases:
::
    prompt%> cd bld
    prompt%> ../src/configure --enable-optimize \
              --enable-shared --disable-debug \
              --with-cubit=/path/to/cubit  \
              --prefix=$HOME/dagmc_bld/CGM
    prompt%> make
    prompt%> make install

HDF5
====

HDF5:  Follow these steps
::
    prompt%> mkdir $HOME/dagmc_bld/HDF5/bld
    prompt%> cd $HOME/dagmc_bld/HDF5
    prompt%> tar xzf ~/hdf5-1.8.7.tar.gz
    prompt%> ln -s hdf5-1.8.7 src
    prompt%> cd bld
    prompt%> ../src/configure --enable-shared --prefix=$HOME/dagmc_bld/HDF5
    prompt%> make
    prompt%> make install


MOAB
````

MOAB: Create a directory to install in
::
    prompt%> mkdir -p $HOME/dagmc_bld/MOAB/bld
    prompt%> cd $HOME/dagmc_bld/MOAB


If installing MOAB from SVN repository:
::
    prompt%> svn co https://svn.mcs.anl.gov/repos/ITAPS/MOAB/trunk
    prompt%> cd trunk
    prompt%> autoreconf -fi
    prompt%> cd ..
    prompt%> ln -s trunk src


If installing MOAB from a tarball, ``MOAB-4.60.tar.gz``:
::
    prompt%> tar xzf ~/MOAB-4.60.tar.gz
    prompt%> ln -s MOAB-4.60 src


In all MOAB cases:
::
    prompt%> cd bld
    prompt%> ../src/configure --enable-optimize \
              --enable-shared --disable-debug \
              --with-cgm=$HOME/dagmc_bld/CGM  \
              --with-hdf5=$HOME/dagmc_bld/HDF5 \
              --prefix=$HOME/dagmc_bld/MOAB
    prompt%> make
    prompt%> make install

**Build and Install MOAB using a Build Script**

The script *build_dagmc_stack.bash* is in UW-Madison's svalinn/DAGMC repository,
`DAGMC git repo <https://github.com/svalinn/DAGMC>`_ under DAGMC/MCNP5/
To see the available options, in the directory containing *build_dagmc_stack.bash* type 
:: 
    prompt%> ./build_dagmc_stack.bash --help

This script is under development.  It is recommended that the components be compiled individually 
at this time.

Applying DAGMC to Specific Monte Carlo Codes
+++++++++++++++++++++++++++++++++++++++++++++

DAG-MCNP5 Build/Install
~~~~~~~~~~~~~~~~~~~~~~~~

If you would like to use DAGMC with MCNP5, known as DAG-MCNP5, you will also need:

* MCNP5.1.51 source code from `RSICC <http://rsicc.ornl.gov>`_
* a local copy of UW-Madison's `DAGMC git repo <https://github.com/svalinn/DAGMC>`_ 

Automatic Installation
=======================

A package has been prepared that includes many of the requires
software libraries and an automated build script.  Because the DAGMC
team is not authorized to distribute `CUBIT
<http://cubit.sandia.gov>`_ nor `MCNP5.1.51 source code
<http://mcnp.lanl.gov>`_, you must acquire those through the
appropriate channels on your own.

Once you have both of those things, you should be able to use the
DagmcBuildPackage to create a working install of DAG-MCNP5.1.51.

Manual Installation
=====================

The following steps are required to install DAG-MCNP5.  Most of these steps are described in more detail below.

1. Install the DAGMC Toolkit as described above
2. Clone a copy of the DAGMC git repo.
3. Apply the appropriate patch from DAGMC/MCNP5/patch/ to your copy of the MCNP5 source code
4. Build & install the patched version of MCNP5

Some assumptions/conventions:

* all operations are carried out in the a subdirectory ``dagmc_bld`` of a user's home directory
* path to CUBIT files is known, e.g. ``/path/to/cubit``
* all tarballs reside in user's home directory
* MCNP5 source code is available in location ``$HOME/dagmc_bld/MCNP5``
* A cloned DAGMC git repo can be found at ``/path/to/DAGMC``

Apply Patch
############

*Apply DAGMC Patch to MCNP5 v1.60*

Perform the following steps:
::
    prompt%> cd $HOME/dagmc_bld/MCNP5/Source
    prompt%> patch -p1 < /path/to/DAGMC/MCNP5/patch/dagmc.patch.5.1.60


Build DAG-MCNP5
################

*Build DAG-MCNP5 from modified code*

One of the easiest ways to build DAG-MCNP5 is directly using the
``makefile`` from the command-line.  To do this, you must know the
``makefile`` options to build a version of MCNP5 without DAGMC,
usually in the form:
::
    prompt%> make build CONFIG="seq plot gfortran" FC=gfortran MARCH=M64``

or similar.  Starting from these options, you can build DAG-MCNP5 from
a patched source code with:
::
    prompt%> make build CONFIG="seq plot gfortran dagmc" FC=gfortran MARCH=M64 \
                 MOAB_DIR=$HOME/dagmc_bld/MOAB CUBIT_DIR=/path/to/cubit/bin \
		 DAGMC_DIR=/path/to/DAGMC/MCNP5/dagmc


If you are less familiar with building MCNP5 from the ``makefile`` you
may want to use the interactive ``install`` script provided by LANL:
::
    prompt%> ./install

Within the ``install`` program you will need to set the DAGMC build options:

* turn on DAGMC mode
* provide the path to MOAB: ``$HOME/dagmc_bld/MOAB``
* provide the path to CUBIT: ``/path/to/cubit``
Your executable should be available as ``$HOME/dagmc_bld/MCNP5/Source/src/mcnp5``.

DAG-Tripoli4 Access
~~~~~~~~~~~~~~~~~~~

Tripoli4 is distributed by CEA/Saclay as a binary executable.  For
access to DAG-Tripoli4, please contact `Jean-Christophe Trama
<mailto:jean-christophe.trama@cea.fr>`_.

FluDAG Build
~~~~~~~~~~~~

FluDAG uses `FLUKA <http://www.fluka.org>`_ from CERN/INFN with the DAGMC Toolkit.
The steps to build and install FluDAG follow.

*Install  FLUKA*

In order to download FLUKA you need to become a registered user, which you can do at 
the FLUKA register page from a link on the main FLUKA page.
Save the user id and password for future FLUKA updates.  
For a 64 bit linux installation the download filename is of the form *fluka20xx.xx-linux-gfor64bitAA.tar.gz*.  
See the `FLUKA <http://www.fluka.org>`_ website in the **Install and Run** item of the *Documentation* menu.

Once the environment variables have been set proceed to the FluDAG installation.

*FluDAG, Alpha Release*
 
Get the FluDAG Alpha Release repository by cloning the UW-Madison DAGMC repository from github.com and checking out the alpha_release branch:
::

    prompt%> cd $HOME
    prompt%> git clone https://github.com/svalinn/DAGMC.git
    prompt%> get checkout alpha_release.

A run-script provided by FLUKA has been modified to allow it to have input filenames
longer than 8 characters.  The following steps assume you have the FLUKA environment
variables defined:
::
 
    prompt%> cd $FLUPRO/flutil
    prompt%> cp rfluka rfluka.orig
    prompt%> cd DAGMC/FluDAG/src
    prompt%> cp rfluka $FLUPRO/flutil

To compile FluDAG you must have the Cubit bin directory in your LD_LIBRARY_PATH.  Add the following
line to your login script:
::

    prompt%> export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/path/to/Cubit/bin

You also need to know where the MOAB libraries are.
These instructions illustrate passing a known MOAB library path to cmake via the -D switch:
::

    prompt%> cd $HOME/DAGMC/FluDAG
    prompt%> mkdir bld
    prompt%> cmake ../src/. -DMOAB_HOME=$HOME/dagmc_bld/MOAB
    prompt%> make

Upon successful compilation the directory *bld* will have the *mainfludag* executable in it.
In order to test FluDAG, an environment variable, named 'FLUDAG', with the path to the *bld* directory must be set:
Add this statement to your login script:
::

    export FLUDAG=${HOME}/DAGMC/FluDAG/bld/

The FluDAG tests are in a separate repository, which can be cloned from the UW-Madison website:
::

    prompt%> git clone https://github.com/svalinn/fludag_testing.git

To run all the tests type:
::

    prompt%> cd $HOME/fludag_testing
    prompt%> ./run_test test_input

Some of the tests are slow, so the above command will take some time.  
If you want to run just the fast tests, or just the magnetic tests
:: 

    prompt%> ./run_test test_fast
    prompt%> ./run_test test_magnetic

The slow tests can be run separately:
::

    prompt%> ./run_test test_slow

At this time, not all tests will pass.

