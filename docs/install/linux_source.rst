.. _linux_source:

========================
Source Install - Linux
========================

------------
Dependencies
------------
PyNE has the following dependencies:

   #. `Fortran compiler <https://gcc.gnu.org/wiki/GFortran>`_
   #. `C++ compiler <https://gcc.gnu.org/>`_
   #. `CMake <http://www.cmake.org/>`_ (>= 2.8.5)
   #. `NumPy <http://www.numpy.org/>`_ (>= 1.8.0)
   #. `SciPy <http://www.scipy.org/>`_
   #. `Cython <http://cython.org/>`_ (>= 0.19.1)
   #. `HDF5 <http://www.hdfgroup.org/HDF5/>`_
   #. `PyTables <http://www.pytables.org/>`_
   #. `Python <http://www.python.org/>`_
   #. `LAPACK <http://www.netlib.org/lapack/>`_
   #. `BLAS <http://www.netlib.org/blas/>`_
   #. `Jinja2 <http://jinja.pocoo.org/>`_

Optional Depenendencies:
   #. `MOAB <https://press3.mcs.anl.gov/sigma/moab-library>`_
   #. `pyMOAB <https://press3.mcs.anl.gov/sigma/moab-library>`_
   #. `DAGMC <http://svalinn.github.io/DAGMC/>`_

To run tutorial and examples:
   #. `jupyter <http://jupyter.org/>`_

Most of the dependencies are readily available through package managers.  Once
all the dependencies are installed, PyNE can be installed. Download and unzip
the source (`zip`_, `tar`_) or checkout a verison from the PyNE repository
(`Github`_).  Then run the following commands from the directory above the
unzipped directory::

    cd pyne/
    python setup.py install --user
    scripts/nuc_data_make

The ``setup.py`` command compiles and installs the PyNE source code.
Note that this command must be done in the top PyNE directory.
The ``nuc_data_make`` builds and installs a database of nuclear data.
Unfortunately, this must be done as a second step because most nuclear 
data is under some form of license restriction or export control which 
prevents the developers from distributing it with PyNE.  However, the 
``nuc_data_make`` program (which is installed by ``setup.py``) will
do its best to find relevant nuclear data elsewhere on your machine
or from public sources on the internet.


Most common install flags:
**************************
The list of the possible flags can be retrieve using:

unzipped directory::
   python setup.py --help


#. optional arguments:
   -h, --help            show this help message and exit
   --clean [CLEAN]       removes the build directory before continuing.
   --user [USER]         Installs into ~/.local

#. cmake:  CMake arguments.

   -D VAR                Set environment variable.
   --build-type BT       Set build type via CMAKE_BUILD_TYPE, e.g. Release
                           or Debug.
   --deps-root DEPS_ROOT The path to the directory containing all
                           dependencies
   --fast  (default)                Will try to compile from assembly, if possible.
                           This is faster than compiling from source.
   --slow                Will NOT try to compile from assembly, if possible.
                           This is slower as it must compile from source.

#. make:  Make arguments.

   -j J                  Degree of parallelism for build.

#. other/dependencies:

   --hdf5 HDF5           Path to HDF5 root directory.
   --moab [MOAB]         Path to MOAB root directory.
   --dagmc [DAGMC]       Path to DAGMC root directory.
   --prefix PREFIX       Prefix for install location.
   --build-dir BLD_DIR   where to place the build directory
   --bootstrap           Bootstraps the PyNE installation, including
                              nuc_data_make and possibly decaygen.


.. _zip: https://github.com/pyne/pyne/zipball/0.5.1
.. _tar: https://github.com/pyne/pyne/tarball/0.5.1
.. _GitHub: http://github.com/pyne/pyne
