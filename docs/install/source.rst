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
   #. `CMake <http://www.cmake.org/>`_ (>= 3.18)
   #. `NumPy <http://www.numpy.org/>`_ (>= 1.8.0)
   #. `SciPy <http://www.scipy.org/>`_
   #. `Cython <http://cython.org/>`_ (>= 0.29.21)
   #. `HDF5 <http://www.hdfgroup.org/HDF5/>`_
   #. `PyTables <http://www.pytables.org/>`_
   #. `Python <http://www.python.org/>`_ (>= 3.10)
   #. `LAPACK <http://www.netlib.org/lapack/>`_
   #. `BLAS <http://www.netlib.org/blas/>`_
   #. `Jinja2 <http://jinja.pocoo.org/>`_

Optional Dependencies:
   #. `MOAB <https://press3.mcs.anl.gov/sigma/moab-library>`_
   #. `DAGMC <https://svalinn.github.io/DAGMC/install/index.html>`_
   #. `OpenMC <https://docs.openmc.org/en/stable/quickinstall.html>`_

To run tutorials and examples:
   #. `jupyter <http://jupyter.org/>`_

Many of these dependencies can be installed using a package manager (i.e., :code:`apt-get` on Linux, `MacPorts <https://www.macports.org/>`_ or `Homebrew <https://brew.sh/>`_ on macOS, or `Anaconda <https://www.anaconda.com/>`_ on both).

---------------------------
Building from Source
---------------------------

With the recent updates to PyNE, the installation process has been streamlined using `pip`_ and `scikit-build-core`_. You can now build PyNE from source by following the steps below:

1. **Download the Source Code**

   Download and unzip the `source code`_, or clone the PyNE repository from `GitHub`_.

2. **Install the Dependencies**

   Make sure that all the build dependencies (Python, C++ compiler, and Fortran compiler) are installed. Many of them can be installed using a package manager.

3. **Install PyNE**

   Run the following command to build and install PyNE:

   .. code-block:: bash

      python -m pip install .

   It is recommended to create a `virtual environment`_ before installing PyNE.

4. **Configure and Build Nuclear Data**

   After installation, you will need to build and install the nuclear data:

   .. code-block:: bash

      nuc_data_make

   This command creates a database of nuclear data. Since much of the data is subject to licensing or export controls, PyNE cannot distribute it directly. Instead, the `nuc_data_make` program attempts to retrieve it from public sources or locate it on your machine.

-----------------------------
Configuring Build Options
-----------------------------

PyNE provides several build options that can be customized based on your needs. These options are set by configuring the `SKBUILD_CMAKE_ARGS` environment variable. Below are the most common options and their default values:

- **DOWNLOAD_HDF5**: Should PyNE download HDF5? (Default: `ON`)
- **ENABLE_EIGEN3**: Should PyNE use Eigen3? (Default: `ON`)
- **DOWNLOAD_EIGEN3**: Should PyNE download Eigen3? (Default: `ON`)
- **ENABLE_LAPACK**: Should PyNE use LAPACK? (Default: `ON`)
- **DOWNLOAD_LAPACK**: Should PyNE download LAPACK? (Default: `ON`)
- **ENABLE_MOAB**: Should PyNE use MOAB? (Default: `ON`)
- **DOWNLOAD_MOAB**: Should PyNE download MOAB? (Default: `ON`)
- **ENABLE_DAGMC**: Should PyNE use DAGMC? (Default: `ON`)
- **DOWNLOAD_DAGMC**: Should PyNE download DAGMC? (Default: `ON`)
- **ENABLE_SPATIAL_SOLVERS**: Should PyNE build advanced AHOT spatial solvers? (Default: `ON`)
- **ENABLE_ENSDF_PROCESSING**: Should PyNE build ENSDF processing tools? (Default: `OFF`)
- **ENABLE_FAST_COMPILE**: Should PyNE use fast compilation? (Default: `OFF`)

To modify these options, export the appropriate flags using `SKBUILD_CMAKE_ARGS`. For example, if you want to disable HDF5 downloading and MOAB/DAGMC, you can set the environment variable as follows:

.. code-block:: bash

    export SKBUILD_CMAKE_ARGS="-DDOWNLOAD_HDF5=OFF; \
                               -DHDF5_ROOT=/path/to/hdf5; \
                               -DDOWNLOAD_EIGEN3=OFF; \
                               -DDOWNLOAD_LAPACK=OFF; \
                               -DENABLE_MOAB=OFF; \
                               -DENABLE_DAGMC=OFF;"

Then, proceed with the installation as described above.

.. _pip: https://packaging.python.org/en/latest/guides/tool-recommendations/
.. _scikit-build-core: https://scikit-build-core.readthedocs.io/en/latest/
.. _virtual environment: https://packaging.python.org/en/latest/guides/installing-using-pip-and-virtual-environments/
.. _source code: https://github.com/pyne/pyne/releases/
.. _GitHub: http://github.com/pyne/pyne
