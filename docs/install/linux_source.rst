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
   #. `Python 2.7 <http://www.python.org/>`_
   #. `LAPACK <http://www.netlib.org/lapack/>`_
   #. `BLAS <http://www.netlib.org/blas/>`_

Optional Depenendencies:
   #. `MOAB <http://trac.mcs.anl.gov/projects/ITAPS/wiki/MOAB>`_
   #. `PyTAPS <https://pythonhosted.org/PyTAPS/index.html>`_

Most of the dependencies are readily available through package managers.  Once
all the dependencies are installed, PyNE can be installed. Download and unzip
the source (`zip`_, `tar`_) or checkout a verison from the PyNE repository
(`Github`_).  Then run the following commands from the unzipped directory::

    cd pyne/
    python setup.py install --user
    scripts/nuc_data_make

The ``setup.py`` command compiles and installs the PyNE source code.
The ``nuc_data_make`` builds and installs a database of nuclear data.
Unfortunately, this must be done as a second step because most nuclear 
data is under some form of license restriction or export control which 
prevents the developers from distributing it with PyNE.  However, the 
``nuc_data_make`` program (which is installed by ``setup.py``) will
do its best to find relevant nuclear data elsewhere on your machine
or from public sources on the internet.


.. _zip: https://github.com/pyne/pyne/zipball/0.4
.. _tar: https://github.com/pyne/pyne/tarball/0.4
.. _GitHub: http://github.com/pyne/pyne
