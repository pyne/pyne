Getting & Installing the DAGMC Toolkit
----------------------------------------

DAGMC is a toolkit that provides direct geometry support to any Monte
Carlo radiation transport code.  The primary code used for development
and testing of DAGMC is `MCNP5 <http://mcnp-green.lanl.gov/>`_,
developed by `Los Alamos National Laboratory <http://www.lanl.gov>`_
and distributed by the `Radiation Safety Information Computing Center
<http://rsicc.ornl.gov>`_.  There has also been experience with MCNPX
(LANL) and Tripoli4 (CEA/Saclay).

These instructions describe the basic steps for downloading and
installing the software libraries that provide the DAGMC toolkit for
integration with Monte Carlo codes.  After this, code-specific
instructions will be give for each code. 

DAGMC Dependencies and User Workflows
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

Manual Geometry Preparation
~~~~~~~~~~~~~~~~~~~~~~~~~~~

Because this stage is highly interactive, most users prefer to perform
this stage on their desktop or laptop computer.  The only software
required for this stage is the interactive Cubit software.

Automated Model Pre-processing
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The primary purpose of this stage is to "translate" from the solid
model representation to the faceted representation.  As such, it
requires Cubit, CGM and MOAB, and typically results in dependencies on
shared libraries.

Model Analysis
~~~~~~~~~~~~~~

During analysis, only the faceted model is necessary, requiring only
MOAB and the physics code.  It is also possible to combine this stage
with the previous one, requiring Cubit, CGM, MOAB and the physics
code, and will also result in dependencies on shared libraries.

Summary of Dependencies
~~~~~~~~~~~~~~~~~~~~~~~

* Some physics package, e.g. MCNP5
   * `MOAB/DAGMC <http://trac.mcs.anl.gov/projects/ITAPS/wiki/MOAB>`_
       * `HDF5 <http://www.hdfgroup.org/HDF5/release/obtain5.html>`_
       * `CGM <http://trac.mcs.anl.gov/projects/ITAPS/wiki/CGM>`_ 
            * ACIS v19, or `CUBIT <http://cubit.sandia.gov>`_ v12.2



Installing the DAGMC Toolkit
++++++++++++++++++++++++++++

The following 4 steps are required to install the MOAB library,
including the DAGMC toolkit, for use in Monte Carlo radiation
transport tools.  Following these steps, all the pre-requisites for
the automated processing stage will be available.

1. Install `CUBIT <http://cubit.sandia.gov>`_ v12.2
2. Install `CGM <http://trac.mcs.anl.gov/projects/ITAPS/wiki/CGM>`_, using the --with-cubit options
3. Install `HDF5 <http://www.hdfgroup.org/HDF5/>`_
4. Install `MOAB <http://trac.mcs.anl.gov/projects/ITAPS/wiki/MOAB>`_,
   using the --with-cgm and --with-hdf5 options (--with-netcdf may
   also be useful but not necessary)

Here are some assumptions/conventions that are used in these instructions:

* all operations are carried out in the a subdirectory ``dagmc_bld`` of a user's home directory
* path to CUBIT files is known, e.g. ``/path/to/cubit``
* all tarballs reside in user's home directory

If these do not apply to you, please modify your steps accordingly.

     *(For a shortcut to installing DAG-MCNP5.1.51 you may be able to use the DagmcBuildPackage .)*

Installing CGM
""""""""""""""

Create a directory to build CGM:
::
    prompt%> mkdir -p $HOME/dagmc_bld/CGM/bld
    prompt%> cd $HOME/dagmc_bld/CGM

If installing from SVN repository:
::
    prompt%> svn co https://svn.mcs.anl.gov/repos/ITAPS/cgm/trunk
    prompt%> cd trunk
    prompt%> autoreconf -fi
    prompt%> cd ..
    prompt%> ln -s trunk src

If installing from a tarball, ``CGM-12.2.0.tar.gz``:
::
    prompt%> tar xzf ~/CGM-12.2.0.tar.gz
    prompt%> ln -s CGM-12.2.0 src

In all cases:
::
    prompt%> cd bld
    prompt%> ../src/configure --enable-optimize --disable-debug \
              --with-cubit=/path/to/cubit  \
              --prefix=$HOME/dagmc_bld/CGM
    prompt%> make
    prompt%> make install

Installing HDF5
""""""""""""""""

Follow these steps
::
    prompt%> mkdir $HOME/dagmc_bld/HDF5/bld
    prompt%> cd $HOME/dagmc_bld/HDF5
    prompt%> tar xzf ~/hdf5-1.8.3.tar.gz
    prompt%> ln -s hdf5-1.8.3 src
    prompt%> cd bld
    prompt%> ../src/configure --prefix=$HOME/dagmc_bld/HDF5
    prompt%> make
    prompt%> make install


Installing MOAB
""""""""""""""""

Create a directory to install MOAB:
::
    prompt%> mkdir -p $HOME/dagmc_bld/MOAB/bld
    prompt%> cd $HOME/dagmc_bld/MOAB


If installing from SVN repository:
::
    prompt%> svn co https://svn.mcs.anl.gov/repos/ITAPS/MOAB/trunk
    prompt%> cd trunk
    prompt%> autoreconf -fi
    prompt%> cd ..
    prompt%> ln -s trunk src


If installing from a tarball, ``MOAB-3.99.tar.gz``:
::
    prompt%> tar xzf ~/MOAB-3.99.tar.gz
    prompt%> ln -s MOAB-3.99 src


In all cases:
::
    prompt%> cd bld
    prompt%> ../src/configure --enable-optimize --disable-debug \
              --with-cgm=$HOME/dagmc_bld/CGM  \
              --with-hdf5=$HOME/dagmc_bld/HDF5 \
              --prefix=$HOME/dagmc_bld/MOAB
    prompt%> make
    prompt%> make install



Applying DAGMC to Specific Monte Carlo Codes
++++++++++++++++++++++++++++++++++++++++++++

Installing DAG-MCNP5
""""""""""""""""""""

If you would like to use DAGMC with MCNP5, known as DAG-MCNP5, you will also need:

* MCNP5.1.51 source code from `RSICC <http://rsicc.ornl.gov>`_
* DAG-MCNP5.1.51 patch file from the UW-Madison

Automatic Installation
'''''''''''''''''''''''

A package has been prepared that includes many of the requires
software libraries and an automated build script.  Because the DAGMC
team is not authorized to distribute `CUBIT
<http://cubit.sandia.gov>`_ nor `MCNP5.1.51 source code
<http://mcnp.lanl.gov>`_, you must acquire those through the
appropriate channels on your own.

Once you have both of those things, you should be able to use the
DagmcBuildPackage to create a working install of DAG-MCNP5.1.51.

Manual Installation
'''''''''''''''''

The following steps are required to install DAG-MCNP5.  Most of these steps are described in more detail below.

1. Install the DAGMC Toolkit as described above
2. Download a copy of the patch file for your version of MCNP:
    * `MCNP5 v1.60 <https://raw.github.com/svalinn/DAGMC/master/MCNP5/patch/dagmc.patch.5.1.60>`_
    * `MCNP5 v1.51 <https://raw.github.com/svalinn/DAGMC/master/MCNP5/patch/dagmc.patch.5.1.51>`_
    * `MCNP4 v1.40 <https://raw.github.com/svalinn/DAGMC/master/MCNP5/patch/dagmc.patch.5.1.40>`_
3. Apply the patch your copy of the MCNP5.1.60 (or appropriate
   version) source code
4. Build & install the patched version of MCNP5

Some assumptions/conventions:

* all operations are carried out in the a subdirectory ``dagmc_bld`` of a user's home directory
* path to CUBIT files is known, e.g. ``/path/to/cubit``
* all tarballs reside in user's home directory
* MCNP5 source code is available in location ``$HOME/dagmc_bld/MCNP5``

Apply DAGMC Patch to MCNP5 v1.60
................................

Perform the following steps:
::
    prompt%> cd $HOME/dagmc_bld/MCNP5/Source
    prompt%> patch -p 0 < ~/dagmc_install.patch.5.1.60


Build DAG-MCNP5 from modified code
....................................

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
                 MOAB_DIR=$HOME/dagmc_bld/MOAB CUBIT_DIR=/path/to/cubit/bin


If you are less familiar with building MCNP5 from the ``makefile`` you
may want to use the interactive ``install`` script provided by LANL:
::
    prompt%> ./install

Within the ``install`` program you will need to set the DAGMC build options:

* turn on DAGMC mode
* provide the path to MOAB: ``$HOME/dagmc_bld/MOAB``
* provide the path to CUBIT: ``/path/to/cubit``
From the main menu, choose ``C`` to compile.

Your executable should be available as ``$HOME/dagmc_bld/MCNP5/Source/src/mcnp5``.

Access to DAG-Tripoli4
"""""""""""""""""""""""

Tripoli4 is distributed by CEA/Saclay as a binary executable.  For
access to DAG-Tripoli4, please contact `Jean-Christophe Trama
<mailto:jean-christophe.trama@cea.fr>`_.

