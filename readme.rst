PyNE: The Nuclear Engineering Toolkit
=====================================
The PyNE project aims to provide a common set of tools for nuclear 
science and engineering needs.

If you are interested in the package itself, or would like to help
and contribute, please let us know either on the mailing list 
(https://groups.google.com/forum/#!forum/pyne-dev, 
pyne-dev@googlegroups.com) or `github`_.

Examples, documentation, and more can be found at 
http://pyne.io/, the official PyNE project site.

.. _github: https://github.com/pyne/pyne

.. install-start

.. _install:

============
Installation
============

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
   #. `Jinja2 <http://jinja.pocoo.org/>`_

Optional Depenendencies:
   #. `MOAB <https://press3.mcs.anl.gov/sigma/moab-library>`_
   #. `PyTAPS <https://pythonhosted.org/PyTAPS/index.html>`_

Additionally, building the documentation requires the following:

   #. `Sphinx <http://sphinx-doc.org/>`_
   #. `sphinxcontrib-bibtex <https://pypi.python.org/pypi/sphinxcontrib-bibtex/>`_
   #. `PrettyTable <https://code.google.com/p/prettytable/>`_
   #. `numpydoc <https://pypi.python.org/pypi/numpydoc>`_
   #. `jupyter <http://jupyter.org/>`_

Most of the dependencies are readily available through package managers. 

------
Binary
------
Binary distributions of the latest release for mac and linux (64-bit) 
using the conda package manager can be installed by running the command::

    conda install -c conda-forge pyne

where VERSION should be replaced with the version number to be installed.

A windows 32-bit binary is also available on conda via the same command but
it is highly experimental and likely broken. Conda binaries do not have 
moab/pytaps/mesh support (yet).

.. _install_source:

------
Source
------
Installing PyNE from source is a two-step process.  First, download and 
unzip the source (`zip`_, `tar`_).  Then run the following commands from 
the unzipped directory::

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

^^^^^^^^^^^^^^^^^^^
Supported Platforms
^^^^^^^^^^^^^^^^^^^
PyNE is currently built and tested daily on the following platforms
on Python 2.7

 * Ubuntu 12.04 - x86_64
 * OSX 10.8 - x86_64

PyNE has pre-built binaries for the lastest release (0.4) on windows
mac and linux 

 * Mac (Python 2.7 and Python 3.3) - x86_64
 * Linux (Python 2.7 and Python 3.3) - x86_64
 * Windows (Python 2.7) - x86

PyNE has known issues on the following platforms

 * Windows (64-bit build currently not feasible)
 * 32-bit platforms (all variants) have known problems - see #315

^^^^^^^^^^^^^^^^^^^^^^^^^^
Conda Install Instructions
^^^^^^^^^^^^^^^^^^^^^^^^^^
After installing anaconda or miniconda from 
`the downloads page <https://www.anaconda.com/distribution/#download-section>`_,
in a new terminal run the following conda install command::

    conda install -c conda-forge pyne

where VERSION should be replaced with the version number to be installed.

If you have any issues, please let us know.

^^^^^^^^^^^^^^^^^^^^^^^^^^
Conda Build Instructions
^^^^^^^^^^^^^^^^^^^^^^^^^^
On mac and linux PyNE can be installed via the package manager conda. 
After installing anaconda or miniconda from 
`the downloads page <https://www.anaconda.com/distribution/#download-section>`_ 
add conda's binary directory to your bash profile by adding::

    export PATH=/path/to/anaconda/bin:$PATH

to your .bashrc or .bash_profile. Then in a new shell::

    conda install conda-build jinja2 nose setuptools pytables hdf5 scipy

on linux you may also need to run::

    conda install patchelf

Then dowload the latest conda-recipes `here 
<https://github.com/conda/conda-recipes/archive/master.zip>`_

cd to the conda-recipes directory and run::

    conda build pyne
    conda install $(conda build --output pyne)
    nuc_data_make

^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
Mac OSX Specific Instructions
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
The simplest method of installing PyNE on mac is via macports. Version 0.4 
can be installed using the following commands(assuming you are using python 2.7)::

   sudo port install py27-pyne
   nuc_data_make --fetch-prebuilt False

The latest development version of PyNE can also be installed from source. The 
instructions below outline how it can be installed
using the homebrew http://brew.sh/ package manager.

Before starting install the command line tools from https://developer.apple.com/downloads/
you will need to create an account in order to download them. After installing brew and
the command line tools run the following commands::

    ruby -e "$(curl -fsSL https://raw.github.com/mxcl/homebrew/go/install)"
    brew doctor
    brew tap homebrew/science
    brew install hdf5
    brew install cmake
    brew install python

Add::

    export PATH=/usr/local/bin:$PATH
    export PATH=/usr/local/share/python:$PATH

to ~/.bash_profile, then::

    source ~/.bash_profile
    sudo pip install numpy
    sudo chown -R $(whoami) /usr/local
    brew install gfortran
    pip install scipy
    pip install cython
    pip install numexpr
    pip install tables

download pyne-staging cd to that directory::

    cd Downloads/pyne-staging
    python setup.py install


Once those lines have been added, run the following command before running 
``nuc_data_make``::

    source ~/.bashrc


.. _zip: https://github.com/pyne/pyne/zipball/0.5.1
.. _tar: https://github.com/pyne/pyne/tarball/0.5.1


^^^^^^^^^^^^^^^^^^^^^^^^^
Ubuntu Build Scripts
^^^^^^^^^^^^^^^^^^^^^^^^^

Scripts for installing PyNE and all its dependencies from scratch on Ubuntu
14.04 - 15.04 are found `here
<https://github.com/pyne/install_scripts/>`_.

.. install-end

============
Contributing
============
We highly encourage contributions to PyNE! If you would like to contribute, 
it is as easy as forking the repository on GitHub, making your changes, and 
issuing a pull request. If you have any questions about this process don't 
hesitate to ask the mailing list (https://groups.google.com/forum/#!forum/pyne-dev, 
pyne-dev@googlegroups.com).
