.. _source:

==============
Source Install 
==============

------------
Dependencies
------------
PyNE has the following dependencies:

   #. `Fortran compiler <https://gcc.gnu.org/wiki/GFortran>`_
   #. `C++ compiler <https://gcc.gnu.org/>`_
   #. `CMake <http://www.cmake.org/>`_ (>= 2.8.5)
   #. `NumPy <http://www.numpy.org/>`_ (>= 1.8.0)
   #. `SciPy <http://www.scipy.org/>`_
   #. `Cython <http://cython.org/>`_ (>= 0.29.21)
   #. `HDF5 <http://www.hdfgroup.org/HDF5/>`_
   #. `PyTables <http://www.pytables.org/>`_
   #. `Python <http://www.python.org/>`_
   #. `LAPACK <http://www.netlib.org/lapack/>`_
   #. `BLAS <http://www.netlib.org/blas/>`_
   #. `Jinja2 <http://jinja.pocoo.org/>`_

Optional Depenendencies:
   #. `MOAB <https://press3.mcs.anl.gov/sigma/moab-library>`_
   #. `DAGMC <https://svalinn.github.io/DAGMC/install/index.html>`__
   #. `OpenMC <https://docs.openmc.org/en/stable/quickinstall.html>`_
   
To run tutorial and examples:
   #. `jupyter <http://jupyter.org/>`_

Many, if not all of these dependencies can be installed using a package manager
(i.e. :code:`apt-get` on Linux, `MacPorts <https://www.macports.org/>`__ or `Homebrew
<https://brew.sh/>`__ on MacOS, or `Anaconda <https://www.anaconda.com/>`__ on both, etc.)


Most common install flags:
**************************
The list of the possible flags can be retrieved using:

  python setup.py --help

#. optional arguments:

  :code:`-h`, :code:`--help`: show this help message and exit

  :code:`--clean`: removes the build directory before continuing.

  :code:`--user`: Installs into :code:`~/.local`

#. cmake:  CMake arguments.

  :code:`-D VAR`: Set environment variable.

  :code:`--build-type BT`: Set build type via :code:`CMAKE_BUILD_TYPE`, e.g. "release"  or "debug".

  :code:`--deps-root DEPS_ROOT`: The path to the directory containing all dependencies

  :code:`--fast`: (default) Will try to compile from assembly, if possible. This is faster than compiling from source.

  :code:`--slow`: Will NOT try to compile from assembly, if possible. This is slower as it must compile from source.

#. make:  Make arguments.

  :code:`-j J`: Degree of parallelism for build.

#. other/dependencies:

  :code:`--hdf5 HDF5_PATH`: Path to HDF5 root direcitory.

  :code:`--moab MOAB_PATH`: Path to MOAB root directory.

  :code:`--dagmc DAGMC_PATH`: Path to DAGMC root directory.

  :code:`--prefix INSTALL_PATH`: Prefix for install location.

  :code:`--build-dir BUILD_DIR`: where to place the build directory

  :code:`--bootstrap`: Bootstraps the PyNE installation, including :code:`nuc_data_make` and possibly decaygen.



Most of the dependencies are readily available through package managers.  Once
all the dependencies are installed, PyNE can be installed. Download and unzip
the source (`zip`_, `tar`_) or checkout a verison from the PyNE repository
(`Github`_).  Then run the following commands from the directory above the
unzipped directory::

    cd pyne/
    python setup.py install --user
    cd
    nuc_data_make

The ``setup.py`` command compiles and installs the PyNE source code.
Note that this command must be done in the top PyNE directory.
The ``nuc_data_make`` builds and installs a database of nuclear data.
This must be done as a second step because most nuclear data is under 
some form of license restriction or export control which prevents the 
developers from distributing it with PyNE.  However, the 
``nuc_data_make`` program (which is installed by ``setup.py``) will
do its best to find relevant nuclear data elsewhere on your machine
or from public sources on the internet.


.. _zip: https://github.com/pyne/pyne/zipball/0.7.1
.. _tar: https://github.com/pyne/pyne/tarball/0.7.1
.. _GitHub: http://github.com/pyne/pyne
